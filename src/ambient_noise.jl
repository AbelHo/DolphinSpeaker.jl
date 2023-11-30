#!/usr/bin/env julia
using SignalAnalysis
using Plots
plotlyjs()
using Statistics

include("PAM.jl")


folname = "/Volumes/Untitled/data/aspod/aspod4/test/2021.03.12_GainTest"

fols = readdir(folname)
arr=Float64[]

for (root, dirs, files) in walkdir(folname)
    for file in files
        try
            if file[end-3:end]==".bin" && file[1]!='.'
                filepath = joinpath(root, file)
                println(filepath)
                data=loadDataBin(filepath)
                println((std(data[:,1:4],dims=1)))
                append!(arr, std(data[:,1]))

                psd(data[:,1:4],fs=400000); ylims!(-60, 20); 
                a=title!(file)
                display(a)
            end
        catch e
            println("FAILED!!!!")
            print(e)
        end
    end
end 

plot(arr/arr[1],marker='.')
xlims!(0,length(arr+5))
ylims!(0,17)

# println("end")



# for fname in fols
#     if fname[end-3:end]==".bin" && fname[1]!="."
#         print(fols,fname)
#     end
# end