#~ cluster with ff
clipsm
clips_max = mapblocks(clipsm; dims=3) do x
    # @info size(x)
    chan = mapslices(energy, x; dims=1) |> vec |> argmax
    x[:, chan]
end

ffts, freqss = compute_rfft(clips_max, fs)#; type=:log)
win = 1:100
i_offset = 100_000
fft_this = @view ffts[:,win .+ i_offset];

plot(freqss, fft_this)
plot(freqss, mapslices(norm_max, fft_this; dims=1))

sum_amp = sum(fft_this; dims=1)' 
sum_amp |> Plots.scatter
Plots.histogram(sum_amp; bins=50)

# try splitting into high and low amplitude groups
high_inds = findall(x -> x >= 4e7, sum_amp)
plot(freqss, fft_this[:,high_inds])
plot(freqss, mapslices(norm_max, fft_this[:,high_inds]; dims=1))

low_inds = findall(x -> x < 4e7, sum_amp)
plot(freqss, fft_this[:,low_inds])
plot(freqss, mapslices(norm_max, fft_this[:,low_inds]; dims=1))


# try clustering
using Clustering
function kmean_cluster(R, k; maxiter=100)
    res_km = kmeans(R, k; maxiter=maxiter, display=:iter)
    for i in 1:k
        inds = findall(res_km.assignments .== i)
        # plot(freqss, R[:,inds]; title="Cluster $i") |> display
    end
    indices = [findall(res_km.assignments .== i) for i in 1:k]
    cluster = selectdim.(Ref(R), Ref(2), indices)
    cluster, indices
    # return [selectdim(R, 2, findall(res_km.assignments .== i)) for i in 1:k]
end

# plot results in one plot inside subplots
function plot_clusters(clus)
    p = Vector{Any}(undef, length(clus))
    for (i, cluster) in enumerate(clus)
        @info "size(cluster): $(size(cluster))"
        # if one of the dimensions is 1, squeeze it
        # find dimensions with size 1
        sz = size(cluster)
        dims_to_squeeze = findall(x -> x == 1, sz)
        @debug "Dims to squeeze: $dims_to_squeeze"
        if !isempty(dims_to_squeeze)
            cluster = dropdims(cluster; dims=dims_to_squeeze[1])
        end
        @debug "After squeeze: $(size(cluster))"
        p[i] = plot(cluster; title="Cluster $i")
        # p[i] = plot(freqss, cluster; title="Cluster $i")
    end
    # plot all clusters stacked vertically
    plot(p...; layout=(length(clus), 1)) |> display
end

sum_amp = sum(fft_this; dims=1)
clus, idx = kmean_cluster(sum_amp, 2; maxiter=100)
plot_clusters(clus)
fft_cluster = [fft_this[:, id] for id in idx]
plot_clusters(fft_cluster)
# choose louder cluster and recluster
biggerAmpl_cluster_ind = argmax(sum.(fft_cluster))
signalonoise_idx = idx[biggerAmpl_cluster_ind]
clus2, idx2 = kmean_cluster(fft_cluster[biggerAmpl_cluster_ind], 2; maxiter=500)
plot_clusters(clus2)

clus, idx = kmean_cluster(fft_this, 2; maxiter=100)
plot_clusters(clus)

k = 2
R = fft_this; #mapslices(norm_max, fft_this; dims=1)'
res_km = kmeans(R, k; maxiter=100, display=:iter)
for i in 1:k
    inds = findall(res_km.assignments .== i)
    plot(freqss, fft_this[:,inds]; title="Cluster $i") |> display
    plot(freqss, mapslices(norm_max, fft_this[:,inds]; dims=1); title="Cluster $i - Norm") |> display
end 



#~ try pca

using MultivariateStats
# R = mapslices(norm_max, fft_this; dims=1)'
## adding amplitude as a feature (FAILED)
sum_R = sum(R;dims=1)
# normalize sum_R
sum_R = (sum_R .- median(sum_R)) |> x-> x ./ maximum(x) .* maximum(R)
R_pca = [R; sum_R]
Plots.heatmap(R_pca)
# plot(R_pca')

# after kmeans for amplitude selection
R_pca = R[:, signalonoise_idx]



#~ setup 
win = 100:200 #1:100
i_offset = 68_000 #100_000
fft_this = @view ffts[:,win .+ i_offset];
R_pca = fft_this
plot(R_pca)

# amplitude feature
sum_amp = sum(fft_this; dims=1)
Plots.scatter(sum_amp')
# normal pca step
pca_dim = 5
M = fit(PCA, R_pca; maxoutdim=pca_dim)
T = MultivariateStats.transform(M, R_pca)
if size(T, 1) != pca_dim
    @warn "Could not get desired PCA dimensions, got $(size(T, 1)) instead of $pca_dim"
    pca_dim = size(T, 1)
    @warn "Adjusting pca_dim to $pca_dim"
end
for i in 1:pca_dim
    Plots.histogram(T[i,:]; bins=50, title="PC $i") |> display
end
# Plots.scatter(T'; title="Each Components")
Plots.plot(T'; title="Each Components")
Plots.scatter(T'; title="Each Components")
# Plots.scatter(sum_R')
M.prinvars #|> @info "PCA variances: "
Plots.plot(M.prinvars; seriestype=:bar, title="PCA Principal Components (Eigenvalues)") |> display
# 

pcas = Array{Any}(undef, pca_dim, 2)
for i in 1:pca_dim
    inds_up = findall(x -> x > 0, T[i,:])
    inds_down = findall(x -> x < 0, T[i,:])
    pcas[i,1] = R_pca[:, inds_up]
    pcas[i,2] = R_pca[:, inds_down]

    p1 = plot(freqss, pcas[i,1]; title="PCA $i high")
    p2 = plot(freqss, pcas[i,2]; title="PCA $i low") 
    plot(p1, p2; layout=(2,1)) |> display
end

## reconstruct original signals from PCA components
y = MultivariateStats.predict(M, R_pca)
R_reconstructed = reconstruct(M, y)

plot(R_reconstructed)
plot(R_pca)

# for i in 1:pca_dim
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
dbs = dbscan(T./2, 1.3e6; min_cluster_size=3, min_neighbors=1)
# dbs = dbscan(R_pca, 1.3e6; min_cluster_size=3, min_neighbors=1)
# dbs = dbscan(log10.(R_pca), 1.5; min_cluster_size=3, min_neighbors=1)
# number of clusters
n_clusters = maximum(assignments(dbs))+1
# plot overlay assignments, plot each column in different color
p=plot(); p2=[plot(;title="DBSCAN Cluster $i") for i in 0:n_clusters-1]
for i in 1:size(R_pca,2)
    plot!(p, R_pca[:,i]; color=assignments(dbs)[i]+1, title="DBSCAN Clustering Results") 
    # plot!(p, log10.(R_pca[:,i]); color=assignments(dbs)[i], title="DBSCAN Clustering Results") #marker=(2,2,:auto,stroke(0)), 
    plot!(p2[assignments(dbs)[i]+1], R_pca[:,i]; color=assignments(dbs)[i]+1)#; title="DBSCAN Cluster $(assignments(dbs)[i]+1)")
end
p |> display
ppp2 = plot(p2...; layout=(length(p2), 1), size=(600, 400 + (n_clusters>4 ? 100*(n_clusters-4) : 0)))

Plots.histogram(assignments(dbs); bins=n_clusters, title="DBSCAN cluster distribution")
pp1 = plot(assignments(dbs); color=assignments(dbs) .+1, xlabel="seconds", ylabel="cluster", title="DBSCAN labels")
pp2 = Plots.scatter(sum_amp'; color=assignments(dbs) .+1)
display(ppp2)
# sleep(2)
ppp3=plot(pp1, pp2; layout=(2,1))


#~ try HDBSCAN
HDBSCAN = MLJ.@load HDBSCAN pkg=MLJScikitLearnInterface
model = HDBSCAN()#(min_cluster_size=3, min_samples=3)
mach = machine(model, R_pca) |> fit! ## wrap model and required data in an MLJ machine and fit
labels = predict_mode(mach, R_pca')
# plot overlay assignments, plot each column in different color
p=plot();
for i in 1:size(R,2)
    plot!(R[:,i]; color=labels[i], title="HDBSCAN Clustering Results") 
    # plot!(log10.(R[:,i]); color=assignments(dbs)[i], title="DBSCAN Clustering Results") #marker=(2,2,:auto,stroke(0)), 
end
plot!()

using UMAP
umap_model = UMAP.UMAP(n_neighbors=15, n_components=2, min_dist=0.1)
embedding = UMAP.fit_transform(umap_model, R)
scatter(embedding[:,1], embedding[:,2]; title="UMAP 2D embedding")  
# try clustering in UMAP space
k = 2
res_km = kmeans(embedding', k; maxiter=100, display=:iter)
for i in 1:k
    inds = findall(res_km.assignments .== i)
    plot(freqss, fft_this[:,inds])
    plot(freqss, mapslices(norm_max, fft_this[:,inds]; dims=1))
end
# try t-SNE
using TSne
tsne_result = tsne(R; dims=2, perplexity=30.0)
scatter(tsne_result[:,1], tsne_result[:,2]; title="t-SNE 2D embedding")  
# try clustering in t-SNE space
k = 2
res_km = kmeans(tsne_result', k; maxiter=100, display=:iter)
for i in 1:k
    inds = findall(res_km.assignments .== i)
    plot(freqss, fft_this[:,inds])
    plot(freqss, mapslices(norm_max, fft_this[:,inds]; dims=1))
end


using MLJ
import MultivariateStats
PCA = MLJ.@load PCA pkg=MultivariateStats
X, y = @load_iris ## a table and a vector

model = PCA(maxoutdim=2)
mach = machine(model, X) |> fit!

Xproj = transform(mach, X)


#~ try
using MLJ

X, labels  = make_moons(400, noise=0.09, rng=1) ## synthetic data with 2 clusters; X
y = map(labels) do label
    label == 0 ? "cookie" : "monster"
end;
y = coerce(y, Multiclass);

DBSCAN = MLJ.@load DBSCAN pkg=Clustering
model = DBSCAN(radius=0.13, min_cluster_size=5)
mach = machine(model)

model = HDBSCAN(min_cluster_size=3, min_samples=3)
mach = machine(model)

## compute and output cluster assignments for observations in `X`:
yhat = MLJ.predict(mach, X)

## get DBSCAN point types:
report(mach).point_types
report(mach).nclusters

## compare cluster labels with actual labels:
compare = zip(yhat, y) |> collect;
compare[1:10] ## clusters align with classes

## visualize clusters, noise in red:
points = zip(X.x1, X.x2) |> collect
colors = map(yhat) do i
   i == 0 ? :red :
   i == 1 ? :blue :
   i == 2 ? :green :
   i == 3 ? :yellow :
   :black
end
using Plots
Plots.scatter(points, color=colors)