# using MultivariateStats
# using Clustering: dbscan, assignments
# using NMF

# #~ setup 
# win = 100:200 #1:100
# i_offset = 68_000 #100_000

# this_window = c_in_ct[216]
# res = detect_impulsetrain2_impl((aufname,data_filt,fs), t[this_window]; res_dir=joinpath(res_dir,"ct_anal2"), display=display, fullplot=true, plot_everytimeintervalhistogram=true)

# this_window = win .+ i_offset
# this_window = c_in_ct[63] #good but dolphin out of fov_angle
this_window = c_in_ct[199]

fft_this = @view ffts[:, this_window];
R_pca = fft_this
plot(R_pca)

# dis = x-> x # display
dis = display
# amplitude feature
sum_amp = sum(fft_this; dims=1)
Plots.scatter(sum_amp')
# normal pca step
pca_dim = 5
M = fit(PCA, R_pca; maxoutdim=pca_dim)
T = MultivariateStats.transform(M, R_pca)
pca_dim_failed = false
if size(T, 1) != pca_dim
    @warn "Could not get desired PCA dimensions, got $(size(T, 1)) instead of $pca_dim"
    pca_dim = size(T, 1)
    @warn "Adjusting pca_dim to $pca_dim"
    pca_dim_failed = true
end
for i in 1:pca_dim
    Plots.histogram(T[i,:]; bins=50, title="PC $i") |> dis
end
# Plots.scatter(T'; title="Each Components")
Plots.plot(T'; title="Each Components")
Plots.scatter(T'; title="Each Components")
# Plots.scatter(sum_R')
M.prinvars #|> @info "PCA variances: "
Plots.plot(M.prinvars; seriestype=:bar, title="PCA Principal Components (Eigenvalues)") |> dis
# 

pcas = Array{Any}(undef, pca_dim, 2)
for i in 1:pca_dim
    local p2
    inds_up = findall(x -> x > 0, T[i,:])
    inds_down = findall(x -> x < 0, T[i,:])
    pcas[i,1] = R_pca[:, inds_up]
    pcas[i,2] = R_pca[:, inds_down]

    p1 = plot(freqss, pcas[i,1]; title="PCA $i high")
    p2 = plot(freqss, pcas[i,2]; title="PCA $i low") 
    plot(p1, p2; layout=(2,1)) |> dis
end
# for i in 1:pca_dimd.png", i_offset))
# savefig(ppp3, @sprintf("temp/dbscan/dbscan_labels_amp_%07d.png", i_offset))
# savefig(ppp4, @sprintf("temp/dbscan/dbscan_zcombined_%07d.png", i_offset))
#     inds = findall(x -> x >= quantile(T[:,i], 0.9) || x <= quantile(T[:,i], 0.1), T[:,i])
#     pcas[i] = R[:, inds]

#     plot(freqss, pcas[i]; title="PCA $i extremes") |> display
#     plot(freqss, mapslices(norm_max, pcas[i]; dims=1); title="PCA $i extremes - Norm") |> display
# end

# try clustering in PCA space
# k = 2
# res_km = kmeans(T, k; maxiter=100, display=:iter)
# for i in 1:k
#     inds = findall(res_km.assignments .== i)
#     plot(freqss, R[:,inds]; title="Cluster $i") |> display
#     plot(freqss, mapslices(norm_max, R[:,inds]; dims=1); title="Cluster $i - Norm") |> display
# end


#~ try KNN
# KNNClassifier = MLJ.@load KNNClassifier pkg=NearestNeighborModels
# # NearestNeighborModels.list_kernels()
# model = KNNClassifier(weights = NearestNeighborModels.Inverse())
# mach = machine(model, X, y) |> fit! ## wrap model and required data in an MLJ machine and fit
# y_hat = predict(mach, X)
# labels = predict_mode(mach, X)

#~ dbscan
dbs = dbscan(T./2, 1.3e6; min_cluster_size=3, min_neighbors=2)
# dbs = dbscan(R_pca, 1.3e6; min_cluster_size=3, min_neighbors=2)
# dbs = dbscan(log10.(R_pca), 1.5; min_cluster_size=3, min_neighbors=1)
# number of clusters
n_clusters = maximum(assignments(dbs))+1
# plot overlay assignments, plot each column in different color
p=plot(); p2=[plot(;title="DBSCAN Cluster $i") for i in 0:n_clusters-1]
for i in 1:size(R_pca,2)
    plot!(p, freqss, R_pca[:,i]; color=assignments(dbs)[i]+1, title="DBSCAN Clustering Results") 
    # plot!(p, log10.(R_pca[:,i]); color=assignments(dbs)[i], title="DBSCAN Clustering Results") #marker=(2,2,:auto,stroke(0)), 
    plot!(p2[assignments(dbs)[i]+1], freqss, R_pca[:,i]; color=assignments(dbs)[i]+1)#; title="DBSCAN Cluster $(assignments(dbs)[i]+1)")
end
p
ppp2 = plot(p2...; layout=(length(p2), 1), size=(600, 400 + (n_clusters>4 ? 100*(n_clusters-4) : 0)),
            legend=false)

Plots.histogram(assignments(dbs); bins=n_clusters, title="DBSCAN cluster distribution")
pp1 = Plots.scatter(t_inS[this_window], assignments(dbs); color=assignments(dbs) .+1, xlabel="seconds", ylabel="cluster", title="DBSCAN labels", markerstrokewidth=0)
pp2 = Plots.scatter(t_inS[this_window], sum_amp'; color=assignments(dbs) .+ 1, markerstrokewidth=0)
ppp3=plot(pp1, pp2; layout=(2,1))
dis(ppp2)
# sleep(2)
ppp4 = plot(ppp2, ppp3; layout=(1,2), size=(1500, 1000))#, title=(pca_dim_failed ? "_failed" : ""))
dis(ppp4)
# savefig(ppp2, @sprintf("temp/dbscan_ct/dbscan_clusters_%07d.png", this_window[1]))
# savefig(ppp3, @sprintf("temp/dbscan_ct/dbscan_labels_amp_%07d.png", this_window[1]))
# savefig(ppp4, @sprintf("temp/dbscan_ct2_minN2_FFT-only2/dbscan_zcombined_ct%05d_c%07d.png", ct_ind, this_window[1]))
# display(ppp4)
pca_dim_failed = false


#~ try NMF
pca_dim = 4
r = nnmf(R_pca, pca_dim)
T = r.H
for i in 1:pca_dim
    Plots.histogram(T[i,:]; bins=50, title="NMF $i") |> dis
end
Plots.plot(T'; title="Each NMF Components")
Plots.scatter(T'; title="Each NMF Components", markerstrokewidth=0)
# try clustering in NMF space
dbs = dbscan(T, 1e3; min_cluster_size=3, min_neighbors=2)
labels = assignments(dbs) .+ 1
# number of clusters
n_clusters = maximum(labels)+1

labels = hdbscan(T', 3) .- 1
n_clusters = labels |> unique |> length
# plot overlay assignments, plot each column in different color
p=plot(); p2=[plot(;title="NMF DBSCAN Cluster $i") for i in 0:n_clusters-1]
for i in 1:size(R_pca,2)
    plot!(p, freqss, R_pca[:,i]; color=labels[i], title="NMF DBSCAN Clustering Results") 
    # plot!(p, log10.(R_pca[:,i]); color=labels[i], title="DBSCAN Clustering Results") #marker=(2,2,:auto,stroke(0)), 
    plot!(p2[labels[i]+1], freqss, R_pca[:,i]; color=labels[i])#; title="DBSCAN Cluster $(assignments(dbs)[i]+1)")
end
p
ppp2 = plot(p2...; layout=(length(p2), 1), size=(600, 400 + (n_clusters>4 ? 100*(n_clusters-4) : 0)),
            legend=false)
Plots.histogram(labels; bins=n_clusters, title="NMF DBSCAN cluster distribution")
pp1 = Plots.scatter(t_inS[this_window], labels; color=labels, xlabel="seconds", ylabel="cluster", title="NMF DBSCAN labels", markerstrokewidth=0)
pp2 = Plots.scatter(t_inS[this_window], sum_amp'; color=labels, markerstrokewidth=0)
ppp3=plot(pp1, pp2; layout=(2,1))
dis(ppp2)
# sleep(2)
ppp4 = plot(ppp2, ppp3; layout=(1,2), size=(1500, 1000))#, title=(pca_dim_failed ? "_failed" : ""))
dis(ppp4)
# savefig(ppp2, @sprintf("temp/dbscan_ct/dbscan_clusters_%07d.png", this_window[1]))
# savefig(ppp3, @sprintf("temp/dbscan_ct/dbscan_labels_amp_%07d.png", this_window[1]))
# savefig(ppp4, @sprintf("temp/dbscan_ct2_minN2_FFT-only2/dbscan_zNMFcombined_ct%05d_c%07d.png", ct_ind, this_window[1]))
# display(ppp4)

## plot ffts vs time
n = size(R_pca, 1)
colors = cgrad(:rainbow, n; rev=true).colors
p = plot()
for i in 1:size(R_pca,1)
    Plots.scatter!(p, R_pca[i,:]; color=colors[i], label=nothing, linewidth=1, alpha=(R_pca[i,:] ./ maximum(R_pca)), markerstrokewidth=0)
end
p

#~ include time into feature
# R = [R_pca; (t_inS[this_window] .- t_inS[this_window[1]])' .* 2e7]
R = [R_pca; (t_inS[this_window] .- t_inS[this_window[1]])' .* (maximum(R_pca)/(t_inS[this_window[end]]-t_inS[this_window[1]]))]#2e7]
plot(R; label=nothing)

# labels = hdbscan(R_pca', 3) .- 1
labels = hdbscan(R', 3) .- 1
n_clusters = labels |> unique |> length
# plot overlay assignments, plot each column in different color
p=plot(); p2=[plot(;title="FFT+time DBSCAN Cluster $i") for i in 0:n_clusters-1]
for i in 1:size(R_pca,2)
    plot!(p, freqss, R_pca[:,i]; color=labels[i], title="NMF DBSCAN Clustering Results") 
    # plot!(p, log10.(R_pca[:,i]); color=labels[i], title="DBSCAN Clustering Results") #marker=(2,2,:auto,stroke(0)), 
    plot!(p2[labels[i]+1], freqss, R_pca[:,i]; color=labels[i])#; title="DBSCAN Cluster $(assignments(dbs)[i]+1)")
end
p
ppp2 = plot(p2...; layout=(length(p2), 1), size=(600, 400 + (n_clusters>4 ? 100*(n_clusters-4) : 0)),
            legend=false)
Plots.histogram(labels; bins=n_clusters, title="NMF DBSCAN cluster distribution")
pp1 = Plots.scatter(t_inS[this_window], labels; color=labels, xlabel="seconds", ylabel="cluster", title="NMF DBSCAN labels", markerstrokewidth=0)
pp2 = Plots.scatter(t_inS[this_window], sum_amp'; color=labels, markerstrokewidth=0)
ppp3=plot(pp1, pp2; layout=(2,1))
dis(ppp2)
# sleep(2)
ppp4 = plot(ppp2, ppp3; layout=(1,2), size=(1500, 1000))#, title=(pca_dim_failed ? "_failed" : ""))
dis(ppp4)
savefig(ppp4, @sprintf("temp/hdbscan_ct2_minClus2_FFTnTime/hdbscan_combined_ct%05d_c%07d.png", ct_ind, this_window[1]))


#~ run all clicks individually
# while true
#     i_offset += 100
#     include("src/dev3_scrap_temp1.jl")
# end
#~ run each click train
# for (ct_idx,now_window) in enumerate(c_in_ct)
#     this_window = now_window
#     ct_ind = ct_idx
#     include("src/dev3_scrap_temp1.jl")
# end
