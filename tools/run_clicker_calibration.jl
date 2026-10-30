# Explicit CLI; including the library itself performs no work.
include(joinpath(@__DIR__,"..","src","clicker_calibration.jl"))

function clicker_cli(args)
    length(args)>=2 || error("Usage: julia --project=. tools/run_clicker_calibration.jl INPUT_ROOT OUTPUT_ROOT [--results=DELIVERY_ROOT] [--observation-policy=strict|terra_estimated_underwater_coordinates] [--device=A|B|both] [--stage=prepare|sync|annotate|fit|review|refine|all] [--annotator=herdr|external] [--notify=telegram|none] [--no-wait] [--no-plots]")
    input,output=args[1:2]
    options=Dict("device"=>"both","stage"=>"all","annotator"=>"herdr","notify"=>"telegram",
        "results"=>"","observation-policy"=>"strict")
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
    options["observation-policy"] in ("strict","terra_estimated_underwater_coordinates") || error("Unsupported observation policy")
    policy=options["observation-policy"]=="strict" ? nothing : Symbol(options["observation-policy"])
    common=(output_dir=output,devices=devices,annotator=Symbol(options["annotator"]),
        notify=Symbol(options["notify"]),wait_annotations=wait_annotations,
        diagnostic_plots=diagnostic_plots,observation_policy=policy)
    result=if isempty(options["results"])
        run_clicker_calibration(input;common...,stage=Symbol(options["stage"]))
    else
        options["stage"]=="all" || error("--results runs the full workflow; use the Julia delivery_only API to deliver existing fits")
        run_full_clicker_calibration(input;common...,results_dir=options["results"])
    end
    println(ClickerCalibration.JSON.json(Dict("status"=>result["status"],"output_dir"=>abspath(output))))
    return result
end

if abspath(PROGRAM_FILE)==@__FILE__
    clicker_cli(ARGS)
end
