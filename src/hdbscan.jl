using MLJ
using Random
using Plots
using PythonCall
using CategoricalArrays
using MLJScikitLearnInterface

## Minimal helper with no try/catch; assumes HDBSCAN is available and exposes `:labels` in fitted_params
## Example usage:
##    X, _labels = make_moons(400, noise=0.09, rng=1)
##    labels = hdbscan(X, 5)
##    Plots.scatter(X.x1, X.x2; color=labels)
function hdbscan(X, min_cluster_size=5; kwargs...)
    # X, _labels = make_moons(400, noise=0.09, rng=1)
    # load HDBSCAN type (assumes available)
    hdb_type = MLJ.@load HDBSCAN pkg=MLJScikitLearnInterface verbosity=0
    model_inst = hdb_type(min_cluster_size=min_cluster_size, kwargs...)
    mach = machine(model_inst, X)
    Base.invokelatest(fit!, mach)# run_hdbscan_minimal()
    fp = fitted_params(mach)

    # Simple extraction: check common places for labels
    labels_raw = nothing
    if haskey(fp, :labels)
        labels_raw = fp[:labels]
    end

    return labels_raw .|> levelcode
end





function run_hdbscan_test()
    ## generate synthetic two-moons data
    X, labels = make_moons(400, noise=0.09, rng=1) ## synthetic data with 2 clusters; X
    y = map(labels) do label
        label == 0 ? "cookie" : "monster"
    end

    model_instance = nothing
    use_model = :none

    ## Try loading HDBSCAN (via MLJScikitLearnInterface). If unavailable, fall back to DBSCAN.
    try
        hdb_type = MLJ.@load HDBSCAN pkg=MLJScikitLearnInterface verbosity=0
        if hdb_type === nothing
            throw(ErrorException("MLJ.@load returned nothing for HDBSCAN"))
        end
        # Use eval on an expression that calls the returned type to avoid world-age issues
        model_instance = eval(:(( $hdb_type )(min_cluster_size=5)))
        use_model = :hdbscan
        println("Using HDBSCAN model (MLJScikitLearnInterface)")
    catch e
        @warn "Could not load HDBSCAN via MLJScikitLearnInterface: $e. Falling back to DBSCAN (Clustering.jl)."
        dbs_type = MLJ.@load DBSCAN pkg=Clustering verbosity=0
        if dbs_type === nothing
            throw(ErrorException("MLJ.@load returned nothing for DBSCAN"))
        end
        model_instance = eval(:(( $dbs_type )(radius=0.13, min_cluster_size=5)))
        use_model = :dbscan
        println("Using DBSCAN model (Clustering)")
    end

    ## Build and fit MLJ machine for the chosen model
    mach = machine(model_instance, X)
    # Use invokelatest to avoid world-age dispatch problems for newly-loaded methods
    Base.invokelatest(fit!, mach)

    ## compute and output cluster assignments for observations in `X`:
    raw_labels = nothing
    if use_model == :hdbscan
        # MLJScikitLearnInterface.HDBSCAN may not implement `predict` through MLJ.
        # The wrapper stores labels in `fitted_params(mach)[:labels]`.
        fp = fitted_params(mach)
        if haskey(fp, :labels)
            raw_labels = fp[:labels]
        else
            # Try to extract an underlying Python object's labels_ attr
            for (_k,v) in pairs(fp)
                try
                    raw_labels = getproperty(v, :labels_)
                    break
                catch
                end
                try
                    raw_labels = PythonCall.getattr(v, "labels_")
                    break
                catch
                end
            end
            # Fallback: find a vector-valued entry with same length as input
            if raw_labels === nothing
                n = length(X.x1)
                for (_k,v) in pairs(fp)
                    try
                        if isa(v, AbstractVector) && length(v) == n
                            raw_labels = v
                            break
                        end
                    catch
                    end
                end
            end
        end
        if raw_labels === nothing
            error("Unable to obtain cluster labels from HDBSCAN fitted parameters; inspect fitted_params(mach).")
        end
    else
        # For DBSCAN via Clustering.jl wrapper, MLJ.predict should work
        raw_labels = MLJ.predict(mach, X)
    end

    ## The predicted labels for clustering models can be Ints, categorical or strings.
    ## Coerce them to a numeric label vector where possible, otherwise keep as-is.
    function labels_to_ints(v)
        # If already numeric, return Int vector
        try
            return Int.(v)
        catch
        end
        # Try parsing after stringifying (works for categorical values like "1", "-1")
        try
            return parse.(Int, string.(v))
        catch
        end
        # Last resort: return a dense mapping from unique values to 1..K
        vals = collect(v)
        uniqs = unique(vals)
        dict = Dict(u => i for (i,u) in enumerate(uniqs))
        return [dict[x] for x in vals]
    end

    clusters = labels_to_ints(raw_labels)

    println("Report summary:")
    try
        println("point_types: ", report(mach).point_types)
        println("nclusters: ", report(mach).nclusters)
    catch
        # some models may not populate these fields
        println("No point_types/nclusters in report for this model.")
    end

    ## compare cluster labels with actual labels (first 10)
    compare = collect(zip(clusters, y))
    println("First 10 cluster vs class:")
    println(compare[1:10])

    ## visualize clusters; treat common noise label -1 as black
    points = zip(X.x1, X.x2) |> collect
    palette = [:red, :blue, :green, :yellow, :purple, :orange, :cyan, :magenta]
    colors = map(clusters) do i
        # If label is -1 (noise in many HDBSCAN/DBSCAN implementations), color black
        if i == -1
            :black
        else
            # Map positive/zero labels to palette with mod1 to avoid out-of-range
            palette[mod1(i, length(palette))]
        end
    end

    scatter([p[1] for p in points], [p[2] for p in points], color=colors, legend=false, title=string("Clustering result (", use_model, ")"))
end


# Helper to convert various label containers to a Vector{Int}
# function _convert_to_int_labels(labels_raw, n)
#     if labels_raw === nothing
#         error("No labels to convert")
#     end
#     # Python-backed arrays
#     if typeof(labels_raw) <: PythonCall.Py
#         builtins = PythonCall.pyimport("builtins")
#         pylist = builtins.list(labels_raw)
#         return [parse(Int, string(x)) for x in pylist]
#     end
#     # Julia integer vector
#     if isa(labels_raw, AbstractVector) && eltype(labels_raw) <: Integer
#         return Int.(labels_raw)
#     end
#     # Categorical vector
#     if isa(labels_raw, CategoricalArrays.CategoricalVector)
#         return CategoricalArrays.levelcode.(labels_raw)
#     end
#     # Try parsing from string representation
#     try
#         return parse.(Int, string.(labels_raw))
#     catch
#     end
#     # Last resort: map unique values to indices
#     vals = collect(labels_raw)
#     uniqs = unique(vals)
#     dict = Dict(u => i for (i,u) in enumerate(uniqs))
#     return [dict[x] for x in vals]
# end

## Minimal OPTICS function: loads and runs OPTICS via MLJScikitLearnInterface and returns integer labels
function optics(X, min_samples=5; kwargs...)
    optics_type = MLJ.@load OPTICS pkg=MLJScikitLearnInterface verbosity=0
    model_inst = eval(:(($optics_type)(min_samples=$(min_samples))))
    mach = machine(model_inst, X)
    Base.invokelatest(fit!, mach)
    fp = fitted_params(mach)

    # Try common places for labels
    labels_raw = nothing
    if haskey(fp, :labels)
        labels_raw = fp[:labels]
    elseif haskey(fp, :fitresult)
        v = fp[:fitresult]
        try
            labels_raw = getproperty(v, :labels_)
        catch
            try
                labels_raw = PythonCall.getattr(v, "labels_")
            catch
            end
        end
    else
        # fallback: pick any vector with same length as input
        n = size(X,2) #length(X.x1)
        for (_k, v) in pairs(fp)
            try
                if isa(v, AbstractVector) && length(v) == n
                    labels_raw = v
                    break
                end
            catch
            end
        end
    end
    labels_raw .|> levelcode
    # return _convert_to_int_labels(labels_raw, length(X.x1))
    # return labels_raw#, length(X.x1)
end

## Minimal test for OPTICS
function run_optics_minimal()
    X, labels = make_moons(300, noise=0.09, rng=3)
    y = map(labels) do label
        label == 0 ? "cookie" : "monster"
    end
    clusters = optics(X, 5)
    n = length(X.x1)
    if length(clusters) != n
        error("OPTICS returned incorrect labels length: $(length(clusters)) != $n")
    end
    compare = collect(zip(clusters, y))
    println("OPTICS first 10 cluster vs class: ")
    println(compare[1:10])
    return clusters
end

## Run test
# run_hdbscan_test()
# a=run_optics_minimal()