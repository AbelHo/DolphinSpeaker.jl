include("run_example.jl")

#~ burst pulse detection
df_fname = "/media/spin/anas2/data/calf/upload/results/results8_20251128_190134/New device/combined_impulse__audacity.txt"
df_clicks = CSV.read(df_fname, DataFrame; header=false)
ct = df_clicks[!, 1]
diff_ct = diff(ct)

Plots.histogram(diff_ct; bins=0.0:0.001:0.1, xlabel="Inter-click interval (s)", ylabel="Counts", title="Inter-click interval histogram")

# get indices of diff_ct between 0.001s and 0.003s
burst_clicks_idx = findall(x -> x >= 0.001 && x <= 0.003, diff_ct) 
burst_clicks_times = ct[burst_clicks_idx .+ 1] # +1
Plots.scatter(burst_clicks_times; xlabel="Time (s)", ylabel="Burst pulse clicks", title="Burst pulse clicks over time", markerstrokewidth=0, alpha=0.5)
# only more than min_click_len clicks in a burst are considered burst pulses
min_click_len = 5

# get burst pulse clicks grouped by bursts
burst_pulses = Vector{Vector{Float64}}()
burst_pulses_trange = Vector{Tuple{Float64, Float64}}()
i = 1
while i <= length(burst_clicks_idx)
    start_idx = burst_clicks_idx[i]
    start_time = ct[start_idx]
    j = i + 1
    while j <= length(burst_clicks_idx) && burst_clicks_idx[j] == burst_clicks_idx[j - 1] + 1
        j += 1
    end
    end_idx = burst_clicks_idx[j - 1] + 1
    end_time = ct[end_idx]
    num_clicks = j - i + 1
    if num_clicks >= min_click_len
        push!(burst_pulses, ct[start_idx:end_idx])
        push!(burst_pulses_trange, (start_time, end_time))
    end
    i = j
end
@info "Detected $(length(burst_pulses)) burst pulses click"
@info "Number of burst pulses: $(length(burst_pulses_trange))"

# export to audacity label file
open( joinpath(dirname(df_fname), "burst_pulses__audacity.txt"), "w") do io
    for (i, (start_t, end_t)) in enumerate(burst_pulses_trange)
        writedlm(io, [start_t end_t "bp_$(i)"], '\t')
    end
end