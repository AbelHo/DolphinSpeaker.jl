
include("audio.jl")
using Plots; plotlyjs()
data, fs = readAudio("/media/spin/anas2/data/marecet/sample/7003.240514151859_465s.flac")
nfft = nextfastfft(round(Int, nfft_inS*fs))
# nfft = 512
win = hanning(nfft)
# win = nothing
specs = mapslices( x -> spectrogram(x, nfft; fs=fs, window=win), data, dims=1)
specs = mapslices( x -> spectrogram(x, nfft, nfft÷4*3; fs=fs, window=win), data, dims=1)
mag_ft = broadcast( x -> pow2db.(x.power) ,specs)
Plots.heatmap(specs[1].time, specs[1].freq, mag_ft[1]; colorbar=false, title="win: $(win === nothing ? "rectangular" : "hamming"), nfft: $nfft")
Plots.heatmap(specs[1].time, log2.(specs[1].freq), mag_ft[1]; colorbar=false, title="win: $(win === nothing ? "rectangular" : "hamming"), nfft: $nfft")


# 3D surface: time x log2(freq) x magnitude
freqss = specs[1].freq
idx = findall(>(0), freqss)  # drop 0Hz
t = specs[1].time
y = freqss[idx]#
y = log2.(freqss[idx])
z = mag_ft[1][idx, :]
# z = specs[1].power[idx, :]
z[findall(z .< 0)] .= 0.001
z = log10.(z) .* 10

z[findall(z .< 0)] .= 0.001
# z = log10.(z) .* 10
# z .= z ./5
# z = z .- minimum(z, dims=2)
# z[findall(z .< 0)] = 0

surface(t, y, z;
    xlabel = "Time (s)",
    ylabel = "Frequency (Hz) (log2)",
    zlabel = "Amplitude (dB)",
    title = "win: $(win === nothing ? "rectangular" : "hamming"), nfft: $nfft",
    colorbar = true,)
    # yticks = freqss[idx])
    # c = :viridis)
    Plots.heatmap(t, y, z)#; c = :grays)#, colorbar=false)

# Plots.heatmap(specs[1].time, specs[1].freq, specs[1].power; colorbar=false)

# Plots.histogram(z |>vec)
# get Statistics of z
# quantile(z|>vec, [0.25, 0.5, 0.75, 0.9, 0.95, 0.99])
# Plots.vline!(quantile(z|>vec, [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]); label="quantiles")

#~ clustering
include("hdbscan.jl")
# time constant
# tc = 0.1 # seconds
# feat_f = log2.(specs[1].freq[idx] ./ tc)
# feat_t = specs[1].time .* tc

# feat_f = log2.(specs[1].freq[idx])
# feat_t = specs[1].time# .* 10
# # quantile(z|>vec, 0.99) - quantile(z|>vec, 0.01)
# feat_z = z |> vec# ./ (quantile(z|>vec, 0.99) - quantile(z|>vec, 0.01))

# Plots.histogram(feat_z)
# Plots.histogram(feat_z[findall(feat_z .> 0)]; labels="z")
# Plots.histogram!(feat_t|>collect; labels="t")
# Plots.histogram!(feat_f; labels="f")

# Plots.heatmap(feat_t,feat_f, reshape(feat_z, size(z)))

feat_f = specs[1].freq[idx] .|> log2
feat_t = specs[1].time .* 100
feat_z = 10 .* log10.(specs[1].power[idx,:])# |> vec
feat_z[findall(feat_z .< 0)] .= 0.001; feat_z = feat_z .|> log10
# feat_z = feat_z .- mean(feat_z; dims=2)# |> vec #(z |> vec)# ./ 5

# 10 .* log10.(specs[1].power[idx,:]) |> eachrow .|> std |> plot
# feat_t |> diff |> plot
# specs[1].freq[idx] .|> log2 |> plot

# try 2
# feat_f = specs[1].freq[idx] .|> log2
# feat_t = 1:length(specs[1].time)
# feat_z = 10 .* log10.(specs[1].power[idx,:]) #|> vec
# feat_z = feat_z .- mean(feat_z; dims=2) |> vec #(z |> vec)# ./ 5

feat_t = (feat_t .- mean(feat_t)) ./ std(feat_t); feat_f = (feat_f .- mean(feat_f)) ./ std(feat_f); feat_z = (feat_z .- mean(feat_z)) ./ std(feat_z)
feat_z = feat_z .* 10
std(feat_f), std(feat_t), std(feat_z)
mean(feat_f), mean(feat_t), mean(feat_z)

R = hcat( repeat(feat_f, size(feat_t,1)) |> vec,
          repeat(feat_t', size(feat_f,1)) |> vec,
          feat_z |> vec )
Plots.heatmap(feat_t, feat_f, reshape(R[:,3], size(z)))
surface(feat_t, feat_f, reshape(R[:,3], size(z)); size=(1400,800))
# thresholding low amplitude points
soft_idx = findall(<(0.001), @view(R[:,3]))
loud_idx = setdiff(1:size(R,1), soft_idx)
R_filt = copy(R); R_filt[soft_idx, :] .= 0.001
Plots.heatmap(t, y, reshape(R_filt[:,3], size(z)))
# Plots.heatmap(t, y, reshape(R[:,3], size(z)))
surface(feat_t, feat_f, reshape(R_filt[:,3], size(z)); size=(1400,800))

# # R[soft_idx, :] .= 0.0
# surface(R[loud_idx,2], R[loud_idx,1], R[loud_idx,3])
# surface(t, y, reshape(R[:,3], size(z)))
# gaussian filter over R
# using ImageFiltering
# R_gf = copy(R)
# kernel = Kernel.gabor(10,10,0.5, 0,5,10,1)
# R_gf[:,3] = imfilter( reshape(R[:,3], size(z)), kernel) |> vec
# Plots.heatmap(t, y, reshape(R_gf[:,3], size(z)))  

# R = R[101:150,:]
labels, m, hd = hdbscan_o(R; min_cluster_size=5, core_dist_n_jobs=Threads.nthreads());#; labels .-= minimum(labels)
# labels = hdbscan(R, 5)
labels |> unique

# R[loud_idx, 3] = R[loud_idx, 3] .|> log10
labels, m, hd = hdbscan_o(R[loud_idx, :]; min_cluster_size=5)#, metric="manhattan")#; labels .-= minimum(labels)
labels_all = ones(Int, size(R,1)) .* -1
labels_all[loud_idx] = labels
labels = labels_all
labels |> unique
# labels, m, hd = hdbscan_o(R, min_cluster_size=3)#; labels .-= minimum(labels)
Plots.heatmap(reshape(labels, size(z)))
Plots.heatmap(reshape(R[:,3], size(z)))

surface(t, freqss[idx], reshape(labels, size(z)); size=(1400,800))

# labels, m, hd = hdbscan_o([R[loud_idx, 1:2] labels./maximum(labels)]; min_cluster_size=5)

# reshape R back to time-frequency matrix
label_mat = reshape(labels, size(z))
Plots.heatmap(label_mat)

label_mat[findall(label_mat .> 0)] .= 0
Plots.heatmap(label_mat)
# Plots.heatmap(t, y, label_mat;
#             #   yticks = freqss[idx],#(y[1:10:end], round.(Int, 2 .^ y[1:10:end])),
#               xlab = "Time (s)",
#               ylab = "Frequency (Hz)",
#               colorbar = false,
#               title = "HDBSCAN Clustering (win: $(win === nothing ? "rectangular" : "hamming"), nfft: $nfft)",
#               )

Plots.heatmap(reshape(R[:,3], size(z)))

# dbscan
dbs = dbscan(R[loud_idx, :], 5; min_cluster_size=5, min_neighbors=3)
# dbs = dbscan(R_pca, 1.3e6; min_cluster_size=3, min_neighbors=2)
# dbs = dbscan(log10.(R_pca), 1.5; min_cluster_size=3, min_neighbors=1)
# number of clusters
n_clusters = maximum(assignments(dbs))+1

#~ try heirarchical clustering
using Clustering
# Rs = R[:,1:52] + R[:,1:52]'
hc = hclust(R+R', linkage=:single)
plot(hc)

D = rand(10, 10)
D += D'
hc = hclust(D, linkage=:single)


hclust_res = hclust(R; linkage=:ward)
labels = cutree(hclust_res, k=5)
Plots.heatmap(reshape(labels, size(z)))
surface(t, freqss[idx], reshape(labels, size(z)); size=(1400,800))

#~ try KNN
using MLJ
KNNClassifier = MLJ.@load KNNClassifier pkg=NearestNeighborModels
X, y = @load_crabs;
# NearestNeighborModels.list_kernels()
model = KNNClassifier(weights = NearestNeighborModels.Inverse())
mach = machine(model, X, y) |> fit! ## wrap model and required data in an MLJ machine and fit
y_hat = predict(mach, X)
labels = predict_mode(mach, X)

# transform labels back to time-frequency matrix
label_mat = reshape(labels, size(z))
Plots.heatmap(t, y, label_mat;
              yticks = freqss[idx],#(y[1:10:end], round.(Int, 2 .^ y[1:10:end])),
              xlab = "Time (s)",
              ylab = "Frequency (Hz)",
              colorbar = false,
              title = "HDBSCAN Clustering (win: $(win === nothing ? "rectangular" : "hamming"), nfft: $nfft)",
              )

###############################################################
# drop the 0 Hz bin (log axis can't show 0)
freqss = specs[1].freq
idx = findall(>(0), freqss)
t = specs[1].time
z = mag_ft[1][idx, :]

# transformed y axis (log2)
y = log2.(freqss[idx])

# choose sensible linear-frequency ticks and clamp to data range
candidate_ticks = [50,100,200,500,1_000,2_000,5_000,10_000,20_000]
freq_min, freq_max = minimum(freqss[idx]), maximum(freqss[idx])
yticks_vals = filter(x -> x >= freq_min && x <= freq_max, candidate_ticks)

# tick positions in the plotted coordinate (log2) and their labels in Hz
yticks_pos = log2.(yticks_vals)
yticks_labels = string.(yticks_vals)

Plots.heatmap(t, y, z;
              yticks = (yticks_pos, yticks_labels),
              xlab = "Time (s)",
              ylab = "Frequency (Hz)",
              colorbar = false)