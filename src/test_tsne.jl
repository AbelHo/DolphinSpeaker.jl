using TSne, Statistics#, MLDatasets

fs_scale = 0.1 # scaling of original sampling rate

rescale(A; dims=1) = (A .- mean(A, dims=dims)) ./ max.(std(A, dims=dims), eps())

alldata, allabels = MNIST.traindata(Float64);
data = reshape(permutedims(alldata[:, :, 1:2500], (3, 1, 2)),
               2500, size(alldata, 1)*size(alldata, 2));
# Normalize the data, this should be done if there are large scale differences in the dataset
X = rescale(data, dims=1);

Y = tsne(X, 2, 50, 1000, 20.0);
Y = tsne(X, 2, 50, 10_000, 20.0);

using Plots
plotlyjs()
theplot = scatter(Y[:,1], Y[:,2], marker=(2,2,:auto,stroke(0)), color=Int.(allabels[1:size(Y,1)]))
Plots.pdf(theplot, "myplot.pdf")

#~ real data
using CSV, DataFrames
label_fname = "/media/spin/anas2/data/marecet/03062024/03062024/0603_093242_Sc1.Table.1.selections.txt"
# label_fname = "/media/spin/Extreme SSD/data/marecet/dataidata/deployment1_18012024_18042024/SDcard3/8338.240214201545.Table.1.selections.txt"
labels = CSV.read(label_fname, DataFrame)
using CategoricalArrays
categories = labels.Notes |> CategoricalArray
labels.id = labels.Notes |> CategoricalArray .|> levelcode

include("media_info.jl")
labels[:,["Begin Time (s)", "End Time (s)"]] = labels[:,["Begin Time (s)", "End Time (s)"]] .- get_duration("/media/spin/anas2/data/marecet/03062024/03062024/Sc1/0603_093242.WAV")/60
# labels = filter(x->x["Begin File"]==basename(aufname), labels)
[categories|>unique labels.id|>unique]

embedding_fname = "/home/spin/Documents/data/temp/0603_093242_WAV_perch-embed.csv"
embedding_fname = "/home/spin/Documents/data/temp/0603_093242_WAV_48000fs__perch-embed.csv"
embedding_fname = "/media/spin/anas/data_res/dolphin/marecet/temp/aqua/0603_103429.WAV_perch-embed.csv"
embedding_fname = "/media/spin/anas/data_res/dolphin/marecet/temp/aqua/0603_103429__48000fs.WAV_perch-embed.csv"
embedding_fname = "/media/spin/anas/data_res/dolphin/marecet/temp/hydromoth/embedding/hydromoth__all_perch-embed.csv"
embedding_fname = "/media/spin/anas/data_res/dolphin/marecet/temp/soundtrap/embedding/fixed__all_perch-embed.csv"


df = CSV.read("/home/spin/Documents/data/temp/20231128_15.16.48_perch-embed-44100.csv", DataFrame)
df = CSV.read("/home/spin/Documents/data/temp/20231128_15.16.48_perch-embed-50000_2.csv", DataFrame)
df = CSV.read("/home/spin/Documents/data/temp/0603_093242_WAV_perch-embed.csv", DataFrame)
df = CSV.read("/home/spin/Documents/data/temp/0603_093242_WAV_48000fs__perch-embed.csv", DataFrame)

df = CSV.read(embedding_fname, DataFrame)
df.start_time = df.start_time * fs_scale; df.end_time = df.end_time * fs_scale
# shift df to the a singular start time
new_start_time = 0.0
for i in 2:size(df,1)
    if df.start_time[i] == 0.0
        new_start_time = df.end_time[i-1]
    end
    df.start_time[i] = df.start_time[i] + new_start_time
    df.end_time[i] = df.end_time[i] + new_start_time
end

df.label = zeros(Int,size(df,1))

# put labels in the df
i_current = 0
for row in eachrow(labels)
    for i = i_current+1:size(df,1)
        if df.start_time[i] < row["Begin Time (s)"] < df.end_time[i]
            df.label[i] = row.id
            i_current = i
            break
        end
    end
end

category_summary = sort(DataFrame([categories|>unique labels.id|>unique [size(filter( x->x.label==i, df),1) for i in labels.id|>unique]], Symbol.(["label", "id", "number"])), "number"; rev=true)
sets = [(filter( x->x.label==i, df)[:,4:end-1] |> Matrix)' for i in category_summary.id]

sets_mean = [mean(sets[i], dims=2) for i in eachindex(sets)]

X = df[:,4:end-1] |> Matrix
Y = tsne(X, 2, 100, 10_000, 20.0);
@time Y = umap(X', 2)'# n_neighbors=3, min_dist=0.6, spread=2.0)' # EXAMPLE: @time res_jl = umap(mnist_x; n_neighbors=10, min_dist=0.001, n_epochs=200)
theplot = scatter(Y[:,1], Y[:,2], marker=(2,2,:auto,stroke(0)))
theplot = scatter(Y[:,1], Y[:,2], marker=(2,2,:auto,stroke(0)), color=df.label, title="t-SNE (all labels)")
savefig(theplot, joinpath(res_dir,"tsne_all.png")); savefig(theplot, joinpath(res_dir,"tsne_all.html"))
theplot = scatter(Y[:,1], Y[:,2], marker=(2,2,:auto,stroke(0)), color=Int.(df.label .> 0), title="t-SNE (combined labels)")
savefig(theplot, joinpath(res_dir,"tsne_combined.png")); savefig(theplot, joinpath(res_dir,"tsne_combined.html"))


z = Y .< [-100 0]
cluster_indices = findall(z[:,1] .&& z[:,2])
z = Y[:,1] .< -4
cluster_indices = findall(z)

z = Y[:,1] .> [100]
cluster_indices = findall(z)

similar_threshold = 1.85
embedding__whistle_centroid = mean(X[cluster_indices,:], dims=1)
similarities = eachrow(X .- embedding__whistle_centroid) .|> norm
plot(similarities); hline!([similar_threshold])
theplot = scatter(Y[:,1], Y[:,2], marker=(2,2,:auto,stroke(0)), color=Int.(similarities .< similar_threshold) )
# theplot = scatter(Y[:,1], Y[:,2], marker=(2,2,:auto,stroke(0)), color=Int.(similarities .< 18), title="t-SNE (whistle centroid (distance<18)")
savefig(theplot, joinpath(res_dir,"tsne_whistle.png")); savefig(theplot, joinpath(res_dir,"tsne_whistle.html"))

jldsave(joinpath(res_dir,"embedding__whistle_centroid_48000fs.jld2"); embedding__whistle_centroid)
write(joinpath(res_dir,"embedding__whistle_centroid_48000fs.csv"), string(embedding__whistle_centroid)[2:end-1])

df.label_perch = similarities .< similar_threshold
event_time = filter( x->x.label_perch, df)[:, [:start_time, :end_time]] |> Matrix
raven_label(event_time, joinpath(res_dir, "raven_embedding-perch_$(basename(aufname)).txt"))

similarity=[@inbounds sum(X[i,:] .* X[i+1,:]) for i in 1:size(X,1)-1]
differences=[@inbounds norm(X[i,:] - X[i+1,:]) for i in 1:size(X,1)-1]
a=plot(similarity, label="similarity")
b=plot(differences, label="differences")

plot(a,b, layout=(2,1))


#~ longer
using JLD2
embeddings = load("/media/spin/anas/data_res/dolphin/marecet/clustering/hydromoth/embedding/hydromoth__all_perch-embed.hdf5")
embeddings = load("/media/spin/anas/data_res/dolphin/marecet/temp/soundtrap/embedding/fixed__all_perch-embed.hdf5")
embeddings = embeddings["default"]'
Y=tsne(embeddings, 2, 100, 10_000, 20.0)


res_dir = "/media/spin/anas/data_res/dolphin/marecet/clustering"
# load previous results
cluster_type = "kmeans"
dd = load("tsne_result.jld2")
X = dd["X"]; Y=dd["Y"]
# kmeans
num_clusters = 6
R = kmeans(Y',num_clusters)
theplot = scatter(Y[:,1], Y[:,2], marker=(2,2,:auto,stroke(0)), color=assignments(R), title="tSNE kmeans=$num_clusters (all labels)", labels=1:6)
dt=5; raven_label([0:dt:(size(X,1)-1)*dt dt:dt:size(X,1)*dt], joinpath(res_dir,"tsne/raven__tsne_kmeans-6.txt"); label_list=assignments(R).|>string )
savefig(theplot, joinpath(res_dir,"tsne/tsne_kmeans-6.png"))
savefig(theplot, joinpath(res_dir,"tsne/tsne_kmeans-6.html"))

theplot = scatter(assignments(R); xlabel="seconds", ylabel="cluster", title="$(cluster_type) labels vs time")
savefig(theplot, joinpath(res_dir,"tsne/tsne_$(cluster_type)_labels.png"))
savefig(theplot, joinpath(res_dir,"tsne/tsne_$(cluster_type)_labels.html"))

# dbscan
cluster_type = "dbscan"
# mkdir(res_dir, "dbscan")
# dbs = dbscan(Y', 10; min_cluster_size=10, min_neighbors=3)
# dbs = dbscan(Y', 5; min_cluster_size=10, min_neighbors=3) # not too good
# dbs = dbscan(Y', 7; min_cluster_size=10, min_neighbors=3)
# dbs = dbscan(Y', 9; min_cluster_size=10, min_neighbors=3)
# dbs = dbscan(Y', 8; min_cluster_size=20, min_neighbors=3)
# dbs = dbscan(Y', 8; min_cluster_size=100, min_neighbors=3)
dbs = dbscan(Y', 8; min_cluster_size=100, min_neighbors=10)
# assignments(dbs) |> plot
theplot = scatter(Y[:,1], Y[:,2], marker=(2,2,:auto,stroke(0)), color=assignments(dbs), title="tSNE DBSCAN (all labels)")
savefig(theplot, joinpath(res_dir,"tsne/tsne_dbscan.png"))
savefig(theplot, joinpath(res_dir,"tsne/tsne_dbscan.html"))

theplot = scatter(assignments(dbs); xlabel="seconds", ylabel="cluster", title="DBSCAN labels vs time")
savefig(theplot, joinpath(res_dir,"tsne/tsne_dbscan_labels.png"))
savefig(theplot, joinpath(res_dir,"tsne/tsne_dbscan_labels.html"))

dt=5; raven_label([0:dt:(size(X,1)-1)*dt dt:dt:size(X,1)*dt], joinpath(res_dir,"tsne/raven__tsne_$(cluster_type).txt"); label_list=assignments(dbs)|>string )

