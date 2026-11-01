# Shared I/O for ClickerCalibration; including this file performs no writes.
const CALIBRATION_ERROR_LOG = Ref{String}(get(ENV, "CLICKER_ERROR_LOG", ""))

read_json(path::AbstractString) = JSON.parsefile(path)
sha256_file(path::AbstractString) = bytes2hex(open(SHA.sha256, path))

function json_safe(value)
    value isa AbstractFloat && return isfinite(value) ? value : nothing
    value isa Symbol && return string(value)
    value isa AbstractDict && return Dict(string(k) => json_safe(v) for (k,v) in value)
    value isa NamedTuple && return Dict(string(k) => json_safe(v) for (k,v) in pairs(value))
    value isa AbstractArray && return map(json_safe, value)
    value isa Tuple && return map(json_safe, collect(value))
    value === missing && return nothing
    return value
end

function write_json(path::AbstractString, value)
    mkpath(dirname(abspath(path)))
    tmp = path * ".tmp." * string(getpid())
    try
        open(tmp, "w") do io
            write(io, JSON.json(json_safe(value), 2), '\n')
            flush(io)
            ccall(:fsync, Cint, (Cint,), Base.fd(io))
        end
        mv(tmp, path; force=true)
    catch
        isfile(tmp) && rm(tmp)
        rethrow()
    end
    return path
end

function log_failure(operation, exception, bt=catch_backtrace())
    path = CALIBRATION_ERROR_LOG[]
    isempty(path) && (path = get(ENV, "CLICKER_ERROR_LOG", joinpath(pwd(), "docs", "error.log")))
    message = sprint(showerror, exception, bt)
    message = replace(message, r"bot[0-9]+:[A-Za-z0-9_-]+" => "bot[REDACTED]")
    record = Dict("timestamp" => string(Dates.now(Dates.UTC)) * "Z",
        "worker" => get(ENV, "HAGENT_NAME", "clicker-calibration"),
        "pid" => getpid(), "operation" => string(operation),
        "exception_type" => string(typeof(exception)), "exception" => message)
    try
        mkpath(dirname(abspath(path)))
        open(path, "a") do io
            write(io, JSON.json(record), '\n')
            flush(io)
            ccall(:fsync, Cint, (Cint,), Base.fd(io))
        end
    catch log_error
        println(stderr, "Could not append calibration failure log: ", log_error)
        println(stderr, "Original failure: ", message)
    end
    return nothing
end

"""Execute a command; capture=true returns unmodified stdout bytes."""
function run_logged(cmd::Cmd; operation, capture::Bool=false)
    err = IOBuffer()
    try
        if capture
            return read(pipeline(cmd; stderr=err))
        end
        run(pipeline(cmd; stderr=err))
        return nothing
    catch exception
        stderr_text = String(take!(err))
        contextual = ErrorException("$(operation): $(cmd)\n$(sprint(showerror, exception))\n$(stderr_text)")
        log_failure(operation, contextual, catch_backtrace())
        throw(contextual)
    end
end

"""
    exercise_error_logging!(path)

Manual pre-launch check of the durable error log; nothing in the calibration workflow calls it.
It appends an intentional marker record, verifies it persisted, moves it to `error_archive.log`, and
then **rewrites** the whole log at `path` without the marker. Run it only before any worker starts,
since concurrent appends during the rewrite would be lost.
"""
function exercise_error_logging!(path::AbstractString)
    CALIBRATION_ERROR_LOG[] = abspath(path)
    marker = "clicker-calibration-error-control-" * string(getpid())
    try
        error(marker)
    catch exception
        log_failure("intentional error logger control", exception, catch_backtrace())
    end
    occursin(marker, read(path, String)) || error("Error logging control was not persisted")
    # The exercised incident is already resolved; retain the exact record in archive.
    lines = readlines(path; keep=true)
    resolved = filter(line -> occursin(marker, line), lines)
    active = filter(line -> !occursin(marker, line), lines)
    archive = joinpath(dirname(path), "error_archive.log")
    open(archive, "a") do io
        for line in resolved
            record = JSON.parse(line)
            record["resolution"] = "Intentional control; verified durable log contains marker before workload launch."
            write(io, JSON.json(record), '\n')
        end
        flush(io)
        ccall(:fsync, Cint, (Cint,), Base.fd(io))
    end
    # Only call this control before workers start, to avoid rewriting concurrent logs.
    open(path, "w") do io
        foreach(line -> write(io, line), active)
    end
    return Dict("marker" => marker, "passed" => true, "archive" => archive)
end
