# Explicit CLI; including the library itself performs no work.
include(joinpath(@__DIR__,"..","src","clicker_calibration.jl"))

function clicker_cli(args)
    length(args)>=2 || error("Usage: julia --project=. tools/run_clicker_calibration.jl INPUT_ROOT OUTPUT_ROOT [--device=A|B|both] [--stage=prepare|sync|annotate|fit|review|refine|all] [--annotator=herdr|external] [--notify=telegram|none] [--no-wait] [--no-plots]")
    input,output=args[1:2]
    options=Dict("device"=>"both","stage"=>"all","annotator"=>"herdr","notify"=>"telegram")
    wait_annotations=true
    diagnostic_plots=true
    for arg in args[3:end]
        if arg=="--no-wait"
            wait_annotations=false
        elseif arg=="--no-plots"
            diagnostic_plots=false
        else
            startswith(arg,"--") && occursin('=',arg) || error("Unknown argument: $arg")
            key,value=split(arg[3:end],'=';limit=2)
            haskey(options,key) || error("Unknown argument: $key")
            options[key]=value
        end
    end
    devices=options["device"]=="both" ? (:A,:B) : (Symbol(options["device"]),)
    result=run_clicker_calibration(input;output_dir=output,devices=devices,
        stage=Symbol(options["stage"]),annotator=Symbol(options["annotator"]),
        notify=Symbol(options["notify"]),wait_annotations=wait_annotations,diagnostic_plots=diagnostic_plots)
    println(ClickerCalibration.JSON.json(Dict("status"=>result["status"],"output_dir"=>abspath(output))))
    return result
end

if abspath(PROGRAM_FILE)==@__FILE__
    clicker_cli(ARGS)
end
