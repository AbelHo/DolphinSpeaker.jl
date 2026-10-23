# Shared loader for the diagnostics scripts: reads a run_contiguous_folders_multitdoa result folder
# (per-file *__cps*.jld2 + the Impulse5/6 webui CSVs) into one table of in-train clicks.
using JLD2, CSV, DataFrames, Statistics, Printf, Random

const D3_IMSIZE = (1520, 2740)                    # (H, W), set_device__ophk_acoustic_D3_clicker
const D3_FOV = (80.0, 54.46086642417594)          # horizontal, vertical (deg)

# inverse of angle2pxd (linear, so exact also outside the frame)
px2az(px; W=D3_IMSIZE[2], fov=D3_FOV[1]) = (px - W / 2) / (W / 2) * fov / 2
py2el(py; H=D3_IMSIZE[1], fov=D3_FOV[2]) = -(py - H / 2) / (H / 2) * fov / 2
vec3(az, el) = (sind(az) * cosd(el), sind(el), cosd(az) * cosd(el))
angd(a, b) = acosd(clamp(sum(a .* b), -1, 1))

# detection result files in pipeline order (sorted by audio file name)
function result_jld2s(res_folder)
    fs_ = filter(f -> occursin("__cps", f) && endswith(f, ".jld2") && !endswith(f, "_angles.jld2"), readdir(res_folder))
    return joinpath.(res_folder, sort(fs_))
end

"""
    load_clicks(res_folder) -> (df, trains)

One row per in-train click over all files: file index, sample, t (s), ppeak, gap_prev/gap_next (ms, within file),
train id, Impulse5 flag, Impulse6 (SRP) az/el, Impulse5 (GCC) az/el. `trains` = per file the train start indices.
"""
function load_clicks(res_folder)
    rows = DataFrame[]
    fs = 0; trains = Vector{Int}[]; aufnames = String[]
    for (k, f) in enumerate(result_jld2s(res_folder))
        j = jldopen(f); tr = j["res_impulsetrain"]; fs = j["fs"]; push!(aufnames, j["res_impulse"].aufname); close(j)
        p = tr.pind_good; n = length(p)
        gp = [i == 1 ? Inf : (p[i] - p[i-1]) / fs * 1000 for i in 1:n]
        gn = [i == n ? Inf : (p[i+1] - p[i]) / fs * 1000 for i in 1:n]
        train = zeros(Int, n); st = [tr.train_start_ind; n + 1]
        for t in 1:length(tr.train_start_ind); train[st[t]:st[t+1]-1] .= t; end
        push!(rows, DataFrame(file=k, sample=p, t=(p .- 1) ./ fs, ppeak=tr.ppeak, gap_prev=gp, gap_next=gn, train=train))
        push!(trains, tr.train_start_ind)
    end
    df = reduce(vcat, rows)
    N = nrow(df)
    rd(key) = CSV.read(joinpath(res_folder, "detection_pixels__color255_webui__$key.csv"), DataFrame)[1:N, :]
    d5 = rd("Impulse5_gccPhatGuided"); d6 = rd("Impulse6_srpPhat")
    df.flag = d5.shape .== "triangle"
    df.az6 = px2az.(d6.px); df.el6 = py2el.(d6.py)
    df.az5 = px2az.(d5.px); df.el5 = py2el.(d5.py)
    df.frame = d6.frame
    return df, (; trains, fs, aufnames)
end

# write a markdown report: title + the captured text in a code block
function write_report(fname, title, body; notes="")
    open(fname, "w") do io
        println(io, "# ", title, "\n")
        isempty(notes) || println(io, notes, "\n")
        println(io, "```\n", body, "```")
    end
    print(body)
end
