module Schedulers

using DataStructures

export
	Scheduler,
	JobID,
	createjob!,
	deletejob!,
	add_dependency!,
	remove_dependency!,
	setfunction!,
	setresult!,
	schedule!,
	cancel!,
	cancelall!,
	iscanceled,
	wantstorun,
	isactive,
	step!,
	status,
	statusstring,
	jobstatus

const JobID = Int
const Timestamp = Int
const Callback = Union{Function,Nothing}

struct JobStatus # TODO: rename
	changedAt::Threads.Atomic{Timestamp}
	runAt::Timestamp # i.e. the value of changedAt when the job was run
end
mutable struct Job
	changedAt::Threads.Atomic{Timestamp} # Only set by Scheduler thread. Can be read by other threads.
	name::String # Only used in error/log messages. Set once.
	f::Union{Function,Nothing}
	spawn::Bool
	edges::Dict{String,JobID} # name=>fromID
	edgesReverse::Set{Tuple{JobID,String}} # Set{(toID,name)}
	result::Any
	status::Symbol # one of :notstarted,:waiting,:spawned,:running,:done
	statusChangedTime::UInt64 # from time_ns()
	runAt::Timestamp
	waitingFor::Set{JobID}
	callbacks::Vector{Function}
end
Job(changedAt::Int, name::String, f::Union{Function,Nothing}, spawn::Bool, result, status::Symbol, runAt::Timestamp) = Job(Threads.Atomic{Int}(changedAt), name, f, spawn, Dict{String,JobID}(), Set{Tuple{JobID,String}}(), result, status, time_ns(), runAt, Set{JobID}(), [])

Job(name::String, f::Function, changedAt::Int, spawn::Bool) = Job(changedAt, name, f,       spawn, nothing, :notstarted, -1)
Job(name::String, result::Any, changedAt::Int, spawn::Bool) = Job(changedAt, name, nothing, spawn, result,  :done,       changedAt)

struct DetachedJob
	jobID::JobID
	runAt::Timestamp
end

struct Scheduler
	threaded::Bool
	timestamp::Ref{Timestamp}
	jobCounter::Threads.Atomic{JobID}
	jobs::Dict{JobID,Job}
	actions::Channel{Function}
	slowActions::Channel{Function}
	dirtyJobs::Vector{JobID}
	activeJobs::Set{JobID}
	detachedJobs::Dict{DetachedJob, UInt64} # -> time_ns() when job started running
	hasSpawned::Ref{Bool}
end
function Scheduler(;threaded=Threads.nthreads()>1)
	threaded && Threads.nthreads()==1 && @warn "Trying to run threaded Scheduler, but threading is not enabled, please set the environment variable JULIA_NUM_THREADS to the desired number of threads."
	Scheduler(threaded, Ref(0), Threads.Atomic{JobID}(1), Dict{JobID,Job}(), Channel{Function}(Inf), Channel{Function}(Inf), [], Set{JobID}(), Dict{DetachedJob,UInt64}(), Ref(false))
end


# --- exported functions ---

function createjob!(s::Scheduler, result::Any, name::Union{String,Nothing}=nothing; spawn::Bool=true)
	id = newjobid!(s)
	name = name==nothing ? "AnonymousJob#$id" : name
	addaction!(s->_createjob!(s,id,result,name,spawn), s)
	id
end
createjob!(f::Function, s::Scheduler, name::Union{String,Nothing}=nothing; kwargs...) = createjob!(s,f,name; kwargs...)
deletejob!(s::Scheduler, jobID::JobID) = addaction!(s->_deletejob!(s,jobID), s)

add_dependency!(   s::Scheduler, dep::Pair{JobID,JobID}, name::String) = addaction!(s->   add_edge!(s,dep,name), s)
remove_dependency!(s::Scheduler, dep::Pair{JobID,JobID}, name::String) = addaction!(s->remove_edge!(s,dep,name), s)

setfunction!(f::Function, s::Scheduler, jobID::JobID) = addaction!(s->_setfunction!(s,jobID,f), s)
setresult!(s::Scheduler, jobID::JobID, result::Any)   = addaction!(s->_setresult!(s,jobID,result), s)

schedule!(callback::Callback, s::Scheduler, jobID::JobID) = addaction!(s->_schedule!(s,jobID,callback), s)
schedule!(s::Scheduler, jobID::JobID)                     = schedule!(nothing,s,jobID)

cancel!(s::Scheduler, jobID::JobID) = addaction!(s->_cancel!(s,jobID),s)
cancelall!(s::Scheduler) = addaction!(s) do s
	for jobID in keys(s.jobs)
		_cancel!(s, jobID)
	end
end

iscanceled(jobStatus::JobStatus) = jobStatus.changedAt[] > jobStatus.runAt

wantstorun(s::Scheduler) = isready(s.actions) || isready(s.slowActions)
function step!(s::Scheduler)
	while isready(s.actions)
		take!(s.actions)(s)
	end
	s.hasSpawned[] = false # TODO: revise this solution
	while !s.hasSpawned[] && isready(s.slowActions) # at most one slow spawning job per step
		take!(s.slowActions)(s)
	end
	updatetimestamps!(s)
end
isactive(s::Scheduler) = !isempty(s.activeJobs)

function status(s::Scheduler)
	# Running (A[2s], B[3s], ...), Detached Running (...), Spawned (C[0.1s], D [...), Waiting (...)
	now = time_ns()
	running  = Tuple{String,Float64}[]
	detached = Tuple{String,Float64}[]
	spawned  = Tuple{String,Float64}[]
	waiting  = Tuple{String,Float64}[]
	for jobID in s.activeJobs
		job = s.jobs[jobID]
		job.status == :running && push!(running,(job.name,(now-job.statusChangedTime)/1e9))
		job.status == :spawned && push!(spawned,(job.name,(now-job.statusChangedTime)/1e9))
		job.status == :waiting && push!(waiting,(job.name,(now-job.statusChangedTime)/1e9))
	end
	for (detachedJob,jobStartTime) in s.detachedJobs
		jobName = "Unknown"
		haskey(s.jobs, detachedJob.jobID) && (jobName = s.jobs[detachedJob.jobID].name)
		push!(detached, (jobName,(now-jobStartTime)/1e9))
	end
	running, detached, spawned, waiting
end
function statusstring(s::Scheduler)
	running, detached, spawned, waiting = status(s)
	strs = String[]
	isempty(running)  || push!(strs, string("Running (",  join(_durationstring.(running),  ", "), ")"))
	isempty(detached) || push!(strs, string("Detached (", join(_durationstring.(detached), ", "), ")"))
	isempty(spawned)  || push!(strs, string("Spawned (",  join(_durationstring.(spawned),  ", "), ")"))
	isempty(waiting)  || push!(strs, string("Waiting (",  join(_durationstring.(waiting),  ", "), ")"))
	join(strs, ", ")
end

"""
	jobstatus(s::Scheduler, jobID::JobID)

Returns one of :doesntexist, :notstarted, :waiting, :spawned, :running, :done.
"""
jobstatus(s::Scheduler, jobID::JobID) = haskey(s.jobs, jobID) ? s.jobs[jobID].status : :doesntexist


# --- internal functions ---

function _setstatus!(s::Scheduler, jobID::JobID, job::Job, status::Symbol, statusChangedTime::UInt64=time_ns())
	if job.status == :running && status==:notstarted
		@warn "Detaching running job $(job.name)"
		s.detachedJobs[DetachedJob(jobID,job.runAt)] = job.statusChangedTime
	end
	wasActive = job.status in (:waiting,:spawned,:running)
	isActive  =     status in (:waiting,:spawned,:running)
	job.status = status
	job.statusChangedTime = statusChangedTime
	wasActive && !isActive &&  pop!(s.activeJobs, jobID)
	isActive && !wasActive && push!(s.activeJobs, jobID)
	nothing
end
_durationstring(t::Tuple{String,Float64}, digits=1) = string(t[1], '[', round(t[2],digits=digits), "s]")

newjobid!(s::Scheduler) = Threads.atomic_add!(s.jobCounter, 1)
newtimestamp!(s::Scheduler) = s.timestamp[]+=1
addaction!(action::Function, s::Scheduler; slow::Bool=false) = (put!(slow ? s.slowActions : s.actions, action); nothing)

function _createjob!(s::Scheduler, jobID::JobID, x, name::String, spawn::Bool)
	job = Job(name, x, newtimestamp!(s), spawn)
	s.jobs[jobID] = job
	nothing
end

function _deletejob!(s::Scheduler, jobID::JobID)
	@assert haskey(s.jobs, jobID) "Trying to delete nonexisting job with id $(jobID)"
	# remove dependencies
	job = s.jobs[jobID]
	for (name,fromID) in job.edges
		remove_edge!(s, fromID=>jobID, name)
	end
	for (toID,name) in job.edgesReverse
		remove_edge!(s, jobID=>toID, name)
	end
	_setstatus!(s, jobID, job, :notstarted) # ensures removal from active jobs and detaching if needed
	delete!(s.jobs, jobID)
	nothing
end

function _setfunction!(s::Scheduler, jobID::JobID, f::Function)
	@assert haskey(s.jobs, jobID) "Trying to update nonexisting job with id $(jobID)"
	job = s.jobs[jobID]
	job.f = f
	setdirty!(s, jobID)
	nothing
end

function _setresult!(s::Scheduler, jobID::JobID, result::Any)
	@assert haskey(s.jobs, jobID) "Trying to update nonexisting job with id $(jobID)"
	job = s.jobs[jobID]
	@assert job.f == nothing
	result==job.result && return nothing # Nothing to do.
	setdirty!(s, jobID)
	updatetimestamps!(s)
	@assert job.status == :notstarted
	job.result = result
	_setstatus!(s, jobID, job, :done)
	nothing
end

_schedule!(callback::Function, s::Scheduler, jobID::JobID, scheduledAt::Timestamp) = _schedule!(s,jobID,callback,scheduledAt)
function _schedule!(s::Scheduler, jobID::JobID, callback::Callback, scheduledAt::Timestamp = s.timestamp[])
	@assert haskey(s.jobs, jobID) "Trying to schedule nonexisting job with id $jobID"
	updatetimestamps!(s)
	job = s.jobs[jobID]
	changedAt = job.changedAt[]
	scheduledAt < changedAt && return # Scheduling out of date.
	slow = !job.spawn

	if job.runAt < changedAt # the job needs to run at a later timestamp than what has been run before
		@assert job.status == :notstarted


		job.runAt = changedAt
		_setstatus!(s, jobID, job, :waiting)
		empty!(job.callbacks)

		# figure out which jobs we are waiting for to finish
		waitingFor = Set{JobID}() # NB: do not reuse old set as it might be referenced from old dependency callbacks!
		for e in job.edges
			fromID = e[2]
			if s.jobs[fromID].result == nothing
				push!(waitingFor, fromID)
				_schedule!(s, fromID, scheduledAt) do x
					delete!(waitingFor,fromID)
					isempty(waitingFor) && addaction!(s->_spawn!(s,jobID,scheduledAt),s; slow=slow)
				end
			end
		end
		job.waitingFor = waitingFor
		@info "$(job.name) waiting for $(length(job.waitingFor)) jobs to finish."
	end
	@assert job.runAt == changedAt
	callback != nothing && push!(job.callbacks, callback)

	isempty(job.waitingFor) && addaction!(s->_spawn!(s,jobID,scheduledAt),s; slow=slow)

	nothing
end


function _setstatustorunning!(s::Scheduler, jobID::JobID, runAt::Timestamp, statusChangedTime::UInt64)
	haskey(s.jobs, jobID) || return # nothing to do if the job was deleted before it finished
	job = s.jobs[jobID]
	runAt < job.changedAt[] && return # nothing to do
	@assert job.status == :spawned
	_setstatus!(s, jobID, job, :running, statusChangedTime)
	nothing
end


function _spawn!(s::Scheduler, jobID::JobID, scheduledAt::Timestamp)
	haskey(s.jobs, jobID) || return # nothing to do if the job was deleted before it finished
	updatetimestamps!(s)
	job = s.jobs[jobID]
	changedAt = job.changedAt[]
	scheduledAt < changedAt && return # Spawned after dependency was finished, but now out of date.

	if job.status == :done
		job.result != nothing && _finish!(s, jobID, job.runAt, job.result, time_ns())
	elseif job.status == :waiting
		@assert isempty(job.waitingFor)
		@info "Spawning $(job.name) ($jobID@$changedAt)"
		s.hasSpawned[] = true
		inputs = getinput(s,jobID)
		_setstatus!(s, jobID, job, :spawned)
		f = job.f
		if s.threaded && job.spawn
			Threads.@spawn runjob!(s, jobID, job.name, job.changedAt, changedAt, f, inputs)
		else
			runjob!(s, jobID, job.name, job.changedAt, changedAt, f, inputs)
		end
	end
	nothing

end

function _cancel!(s::Scheduler, jobID::JobID)
	@assert haskey(s.jobs, jobID) "Trying to cancel nonexisting job with id $jobID"
	updatetimestamps!(s)
	job = s.jobs[jobID]
	job.result == nothing && setdirty!(s, jobID) # Do not invalidate jobs that have finished!
end

function _finish!(s::Scheduler, jobID::JobID, runAt::Timestamp, result::Any, statusChangedTime::UInt64)
	haskey(s.jobs, jobID) || return # nothing to do if the job was deleted before it finished
	updatetimestamps!(s)
	job = s.jobs[jobID]
	changedAt = job.changedAt[]
	if changedAt == runAt # we finished the last version of the job
		@assert job.status in (:running,:done)
		job.result = result
		_setstatus!(s, jobID, job, :done, statusChangedTime)
		if !(result isa Exception)
			for callback in job.callbacks
				callback(job.result)
			end
		end
		empty!(job.callbacks)
	else
		# remove from list of detached jobs
		pop!(s.detachedJobs, DetachedJob(jobID,runAt))
	end

	result isa Exception && (@warn result; throw(result))

	nothing
end

function getinput(s::Scheduler, jobID::JobID)
	job = s.jobs[jobID]
	Dict{String,Any}((name=>s.jobs[inputJobID].result) for (name,inputJobID) in job.edges)
end

function add_edge!(s::Scheduler, dep::Pair{JobID,JobID}, name::String)
	fromID,toID = dep
	@assert haskey(s.jobs, toID)   "Trying to add edge to nonexisting job with id $(toID)"
	@assert haskey(s.jobs, fromID) "Trying to add edge from nonexisting job with id $(fromID)"
	from,to = s.jobs[fromID], s.jobs[toID]
	@assert !haskey(to.edges, name) # we don't allow adding edges that already exist
	to.edges[name] = fromID
	push!(from.edgesReverse, (toID,name))
	setdirty!(s, toID)
	nothing	
end
function remove_edge!(s::Scheduler, dep::Pair{JobID,JobID}, name::String)
	fromID,toID = dep
	@assert haskey(s.jobs, toID)   "Trying to remove edge to nonexisting job with id $(toID)"
	@assert haskey(s.jobs, fromID) "Trying to remove edge from nonexisting job with id $(fromID)"
	from,to = s.jobs[fromID], s.jobs[toID]
	# we don't allow removing edges that doesn't exist
	@assert haskey(to.edges, name)
	@assert to.edges[name]==fromID
	delete!(to.edges, name)
	delete!(from.edgesReverse, (toID,name))
	setdirty!(s, toID)
	nothing
end

function setdirty!(s::Scheduler, jobID::JobID)
	job = s.jobs[jobID]
	job.changedAt[] = newtimestamp!(s)
	push!(s.dirtyJobs, jobID)
	nothing
end

function updatetimestamprec!(s::Scheduler, jobID::JobID)
	haskey(s.jobs, jobID) || return # happens if the job was deleted
	job = s.jobs[jobID]
	changedAt = job.changedAt[]

	_setstatus!(s, jobID, job, :notstarted)
	job.result = nothing

	for (toID,name) in job.edgesReverse
		toJob = s.jobs[toID]
		Threads.atomic_max!(toJob.changedAt,changedAt) < changedAt && updatetimestamprec!(s,toID)
	end
end

function updatetimestamps!(s::Scheduler)
	while !isempty(s.dirtyJobs) # go through in reverse order to minimize changes to dependency graph
		updatetimestamprec!(s,pop!(s.dirtyJobs))
	end
end

function runjob!(s::Scheduler, jobID::JobID, jobName::String, changedAt::Threads.Atomic{Timestamp}, runAt::Timestamp, f::Function, input::Dict{String,Any})
	jobStatus = JobStatus(changedAt, runAt)
	result = nothing
	if !iscanceled(jobStatus) # ensure the job was not changed after it was scheduled
		try
			@info "Running $jobName ($jobID@$runAt) with $(length(input)) parameter(s) in thread $(Threads.threadid())."
			jobStartTime = time_ns()
			addaction!(s->_setstatustorunning!(s,jobID,runAt,jobStartTime), s)

			# sanity check input and propagate errors
			for (name,v) in input
				@assert v != nothing "Job \"$jobName\" input \"$name\" has not been computed."
				@assert !(v isa Exception) "Job \"$jobName\" input \"$name\" errored: $v."
			end
			result = f(jobStatus, input)
			dur = round((time_ns()-jobStartTime)/1e9,digits=1)
			isCanceled = iscanceled(jobStatus)
			detachedStr = isCanceled ? " (detached job)" : ""
			@info "Finished running$detachedStr $jobName[$(dur)s] ($jobID@$runAt) in thread $(Threads.threadid())."
			isCanceled || @assert result!=nothing "Job $jobName returned nothing"
		catch err
			@warn "Job $jobName ($jobID@$runAt) failed"
			showerror(stdout, err, catch_backtrace())
			result = err
		end
	end
	statusChangedTime = time_ns()
	addaction!(s->_finish!(s,jobID,runAt,result,statusChangedTime), s)
	result isa Exception && throw(result)
	nothing
end

end
