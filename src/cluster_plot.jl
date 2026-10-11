fft_this = @view ffts[:, this_window];
R_pca = fft_this
plot(R_pca)

# clusttering_method = :OPTICS #:HDBSCAN
# clusttering_method = :HDBSCANcondensed_tree_
# dis = x-> x # display
dis = display
# amplitude feature
sum_amp = sum(fft_this; dims=1)

#~ include time into feature
# R = R_pca
# R = [R_pca; (t_inS[this_window] .- t_inS[this_window[1]])' .* 2e7]
# add time scaled to max of R_pca as feature
R = [R_pca; (t_inS[this_window] .- t_inS[this_window[1]])' .* (maximum(R_pca)/(t_inS[this_window[end]]-t_inS[this_window[1]]))]#2e7]
# add total amplitude as feature
# R = [R_pca;
#     (t_inS[this_window] .- t_inS[this_window[1]])' .* (maximum(R_pca)/(t_inS[this_window[end]]-t_inS[this_window[1]])); #2e7]
#     sum_amp ./ (maximum(sum_amp) / maximum(R_pca))
#     ]

# save R as HDF5
# using HDF5
# h5open(joinpath(outdir, @sprintf("%s_ct%05d_c%07d_features.h5", clusttering_method, ct_ind, this_window[1])), "w") do h5
#     write(h5, "R", R)
# end

plot(R; label=nothing)

# labels = hdbscan(R_pca', 3) .- 1
labels = nothing
if clusttering_method == :HDBSCAN
    labels = hdbscan(R', 3); labels .-= minimum(labels)
    # labels, m, hd = hdbscan(R', 3); labels .-= minimum(labels)
elseif clusttering_method == :OPTICS
    labels = optics(R', 3); labels .-= minimum(labels)
elseif clusttering_method == :HDBSCAN_O
    labels, clusterer, df_clustree = hdbscan_o(R'; min_cluster_size=3); labels .-= minimum(labels)
end

n_clusters = labels |> unique |> length
# plot overlay assignments, plot each column in different color
p=plot(); p2=[plot(;title="FFT+time $clusttering_method Cluster $i") for i in 0:n_clusters-1]
for i in 1:size(R_pca,2)
    plot!(p, freqss, R_pca[:,i]; color=labels[i], title="$clusttering_method Clustering Results") 
    # plot!(p, log10.(R_pca[:,i]); color=labels[i], title="DBSCAN Clustering Results") #marker=(2,2,:auto,stroke(0)), 
    plot!(p2[labels[i]+1], freqss, R_pca[:,i]; color=labels[i])#; title="DBSCAN Cluster $(assignments(dbs)[i]+1)")
end
p
ppp2 = plot(p2...; layout=(length(p2), 1), size=(600, 400 + (n_clusters>4 ? 100*(n_clusters-4) : 0)),
            legend=false)#, xtick=:none, title="", titlefont=Plots.font(4))
Plots.histogram(labels; bins=n_clusters, title="$clusttering_method cluster distribution")
pp1 = Plots.scatter(t_inS[this_window], labels; color=labels, xlabel="seconds", ylabel="cluster", title="$clusttering_method labels", markerstrokewidth=0)
pp2 = Plots.scatter(t_inS[this_window], sum_amp'; color=labels, markerstrokewidth=0, alpha=0.4)
rgbs, rgba = sig2rgb(clips_max[:, this_window]; fs=fs, rgb_bands=rgb_bands)
rgbs[4,:] .+= rgbs_alpha_offset
colorings = rgbs|> eachcol .|> x-> RGBA(x...)
pp3 = Plots.scatter(t_inS[this_window], sum_amp'; title="Frequency Color", color=colorings, bg=:black,markerstrokewidth = 0, alpha=rgbs[4,:])#, hover = string.(this_window))

ppp3=plot(pp1, pp2, pp3; layout=(3,1), legend=false)
dis(ppp2)
# sleep(2)
ppp4 = plot(ppp2, ppp3; layout=(1,2), size=(1500, 1000))#, title=(pca_dim_failed ? "_failed" : ""))
dis(ppp4)
savefig(ppp4, joinpath(outdir, @sprintf("%s_combined_ct%05d_c%07d.png", clusttering_method, ct_ind, this_window[1])))





#~ run each click train
# this_window=nothing; ct_ind=nothing
# outdir = "temp/hdbscan_ct2_minClus3_FFTnTime_2"
# mkpath(outdir)
# open( joinpath(outdir, "hdbscan_clusters.csv"), "w") do io
#     writedlm(io, ["c_ind" "time_inS" "ct_ind" "cluster_label"], ',')
# end
# for (ct_idx,now_window) in enumerate(c_in_ct)
#     this_window = now_window
#     ct_ind = ct_idx
#     include("src/cluster_plot.jl")
#     CSV.write(joinpath(outdir, "hdbscan_clusters.csv")),
#         DataFrame(c_ind=this_window, time_inS=t_inS[this_window], ct_ind=ct_ind, cluster_label=labels;
#         append=true)
#     open( joinpath(outdir, "hdbscan_clusters.csv"), "a") do io
#         writedlm(io, [this_window t_inS[this_window] (ones(size(this_window)) .* ct_idx) labels], ',')
#     end
# end

# function plot_scatters(x, y1, y2; kwargs...)
#     pp1 = Plots.scatter(x, y1; color=y1, xlabel="seconds", ylabel="cluster", title="$clusttering_method labels", markerstrokewidth=0)
#     pp2 = Plots.scatter(x, y2; color=y1, markerstrokewidth=0, alpha=0.4)
#     ppp3=plot(pp1, pp2; layout=(2,1), kwargs...)
# end
# i=138
# for i in 1:length(c_in_ct)
#     win = c_in_ct[i]
#     plot_scatters(t_inS[win], labels[win], sum_amp[win]; size=(1000,600), title="Cluster $i") |> display
# end

