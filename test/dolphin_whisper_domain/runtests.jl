# Dedicated original-checkout domain entry point; existing test/runtests.jl is untouched.
# JULIA_DEPOT_PATH=/writable/private/depot:/readonly/package/depot JULIA_NUM_THREADS=1 \
#   julia --startup-file=no --history-file=no --compiled-modules=no --project=. \
#   test/dolphin_whisper_domain/runtests.jl
using Test, SHA, TOML, Dates, JSON
const DOMAIN_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DOMAIN_OUTPUT = haskey(ENV, "DW_DOMAIN_OUTPUT") ? abspath(ENV["DW_DOMAIN_OUTPUT"]) : mktempdir(prefix="dolphin-domain-")
mkpath(DOMAIN_OUTPUT)
const DOMAIN_ERROR_LOG = get(ENV, "DW_DOMAIN_ERROR_LOG", joinpath(DOMAIN_OUTPUT, "error.log"))
function domain_failure(error, trace; operation)
    record = Dict("timestamp"=>string(now(UTC)), "worker"=>"dw-localization-1023a",
        "pid"=>getpid(), "operation"=>operation, "exception_type"=>string(typeof(error)),
        "message"=>sprint(showerror,error), "traceback"=>sprint(showerror,error,trace))
    try
        mkpath(dirname(DOMAIN_ERROR_LOG))
        open(DOMAIN_ERROR_LOG,"a") do io
            println(io,JSON.json(record));flush(io)
        end
    catch logging_error
        @error "Domain error logging failed; preserving original exception" exception=(logging_error,catch_backtrace())
    end
end
function verify_domain_sources()
    specification = TOML.parsefile(joinpath(@__DIR__,"source_hashes.toml"))
    actual = Dict{String,String}()
    for (relative,expected) in specification["files"]
        value=bytes2hex(sha256(read(joinpath(DOMAIN_ROOT,relative))))
        value==expected || error("Pinned source changed: $relative; expected $expected, found $value. Review intentionally; never adopt a newer hash silently.")
        actual[relative]=value
    end
    actual
end
try
    VERSION==v"1.12.6" || error("This domain receipt requires Julia 1.12.6; got $VERSION")
    Threads.nthreads()==1 || error("Run with JULIA_NUM_THREADS=1")
    global DOMAIN_SOURCE_HASHES=verify_domain_sources()
    include(joinpath(@__DIR__,"suite.jl"))
    verify_domain_sources()==DOMAIN_SOURCE_HASHES || error("Source changed during domain tests")
catch error
    domain_failure(error,catch_backtrace();operation="original Julia domain test runner")
    rethrow()
end
