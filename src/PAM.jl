using FileIO 

#load acoustic file
function loadDataBin(fname;channels=5, duration=nothing,fs=400000,datatype=Int16)
    if duration == nothing # load the entire file
        duration = stat(fname).size/channels/fs/sizeof(datatype) #16 bits per sample, 2 bytes
    end
    data = Array{datatype}(undef, Int(floor(fs*duration*channels)))
    read!(fname, data)
    data = transpose( reshape(data,channels,size(data,1)÷channels) )
end

function loadDataBin2(fname; skiplen_inS=0, channels=5, duration=nothing,fs=400000,datatype=Int16)
    f = open(fname)
    skip = skiplen_inS*channels*fs*datatype.size 
    skip = skip-mod(skip,channels*datatype.size) |> Int # round it down to the beginnning of the 1st channel instead of being in between channels
    
    seek(f, skip)

    if duration == nothing # load the entire file
        duration = (stat(fname).size - skip)/channels/fs/datatype.size #16 bits per sample, 2 bytes
    end
    data = Array{datatype}(undef, Int(fs*duration*channels))
    read!(f, data)
    data = transpose( reshape(data,channels,size(data,1)÷channels) )
    close(f)
    return(data)
end

function loadDataBinEndFile(fname;channels=5, duration=nothing,fs=400000,datatype=Int16, arr_type=Array)
    f = open(fname)
    # println(position(f))
    if duration == nothing # load the entire file
        duration = stat(fname).size/channels/fs/datatype.size #16 bits per sample, 2 bytes
    end
    posi = Int(round(stat(f).size-duration*fs*channels*datatype.size));
    seek(f, posi)
    println(position(f))
    # data = datatype.(read(f,Int(fs*duration*channels)))

    data = arr_type{datatype}(undef, Int(fs*duration*channels))
    read!(f, data)
    close(f)
    data = transpose( reshape(data,channels,size(data,1)÷channels) )
    return data, posi
end









function loadDataBinEnd(fol;channels=5, duration=nothing,fs=400000,datatype=Int16, arr_type=Array)
    fols = readdir(fol)
    fname = joinpath(fol,fols[end])
    println(fname)
    
    f = open(fname)
    # println(position(f))
    if duration == nothing # load the entire file
        duration = stat(fname).size/channels/fs/datatype.size #16 bits per sample, 2 bytes
    end
    posi = Int(round(stat(f).size-duration*fs*channels*datatype.size));
    try
        seek(f, posi)        
    catch e
        posi<0 || (posi=0; println("posi=0"); sleep(2) )
    end
    println(position(f))
    # data = datatype.(read(f,Int(fs*duration*channels)))

    data = arr_type{datatype}(undef, Int(fs*duration*channels))
    try
        read!(f, data)
    catch e
        println("************** cant read data, wait.....")
        sleep(1)
        read!(f,data)
    end
    close(f)
    data = transpose( reshape(data,channels,size(data,1)÷channels) )
    return data, posi, fname
end

## for continuous analysis
function loadDataBinEndRegular(fol;channels=5, duration=nothing,fs=400000,datatype=Int16, arr_type=Array)
    fols = readdir(fol)
    fname = joinpath(fol,fols[end])
    println(fname)
    
    f = open(fname)
    # println(position(f))
    if duration == nothing # load the entire file
        duration = stat(fname).size/channels/fs/datatype.size #16 bits per sample, 2 bytes
    end
    winlen=duration*fs*channels*datatype.size
    posi = round(stat(f).size-winlen);
    posi = posi-mod(posi,winlen)
    try
        seek(f, posi)        
    catch e
        posi<0 || (posi=0; println("posi=0"); sleep(2) )
    end
    println(position(f))
    # data = datatype.(read(f,Int(fs*duration*channels)))

    data = arr_type{datatype}(undef, Int(fs*duration*channels))
    try
        read!(f, data)
    catch e
        println("************** cant read data, wait.....")
        sleep(1)
        read!(f,data)
    end
    # close(f)
    data = transpose( reshape(data,channels,size(data,1)÷channels) )
    return data, posi, f
end

function loadDataContinuous(file, data)
    if eof(file)
        #next file
        close(file)
        fol = dirname(file.name[7:end-1])
        fols = readdir(fol)
        fname = joinpath(fol,fols[end])
        file=open(fname)
        println(fname)
    end    
    println(position(file))
    waitcounter=0;
    while stat(file).size < position(file)+4000000
        # println("wait...")
        waitcounter+=1
        sleep(0.1)
    end
    read!(file,data)
    return data,file,waitcounter
end


function loadAcoustic(fname;duration=nothing,volt=false)
    data = loadDataBin(fname,channels=5,duration=duration,fs=400000,datatype=Int16)
    if volt
        data = data*10/32768 # # 16 bits per sample, 1 bit for ± sign 2^15=32768; ±10V measurement
    end
    return data
end

function loadAcousticPressure(fname;duration=nothing)
    data = loadDataBin(fname,channels=5,duration=duration,fs=400000,datatype=Int16)

    ## 16 bits per sample, 1 bit for ± sign 2^15=32768; ±10V measurement
    ## -211dB hydrophone 
    ## 200x default gain 
    ## gain setting 1111 (max): x16V/V
    ## data = data*10/32768 * 10^(211/20) /200 /16
    ## 10/32768 * 10^(211/20) /200 /16 -> 3383.7641642911535
    # return in µPa 
    return data*3383.7641642911535
end

function loadFlac(fname)
    # requires FLAC
    load(fname)
end

########## compass/AHRS

function loadCompass(fname,duration=nothing)
    #1:3 yaw,roll,pitch;  4:6 accelX,Y,Z     7:9 gyroX,Y,Z
    #10  temperature(C)    11 sensor start time  12 epoch time from computer(python)
    label = ["yaw","roll","pitch","accelX","accelY","accelZ","gyroX","gyroY","gyroZ",
    "temperature","time_sensor","time_comp"]
    data = loadDataBin(fname,channels=12,duration=duration,fs=10,datatype=Float64)
    data[:,end] = data[:,end] .+ 28800 # UTC+8 for Singapore local time, python auto converts but julia doesn't. 8*60*60 = 28800
    return data,label
end

function readCSV(fname;tailrows=20,selfheader="date,temperature,pressure,humidity",datefmt="yyyymmdd_HHMMSS")
    tailrows = string(tailrows)
    a=CSV.read(IOBuffer("$selfheader\n"
            *read(`tail -n $tailrows $fname`, String)),
        dateformat=datefmt)
end


################  ANALYSIS ################
using DSP
using SignalBase
function psd2(data; fs=1.0, nfft=512, noverlap=div(nfft,2),
    window=hamming(nfft), xscale=:auto, yrange=50)
    p=[];pow=Array{Float64}(undef, Int(nfft/2+1), size(data,2));
    for i = 1:size(data,2)
        p = welch_pgram(data[:,i], nfft, noverlap; fs=inHz(fs), window=window)
        pow[:,i] = 10*log10.(p.power)
    end

    return pow, p.freq
end
