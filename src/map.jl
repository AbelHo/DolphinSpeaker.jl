using FileIO
using Geodesy: LLA, euclidean_distance, ENUfromLLA, wgs84
using ARLToolkit
import ARLToolkit.geo2pos
import Statistics.mean
# export mean

AOI(LLA(20.734398+.03, -156.903779-.03), LLA(20.734398-.03, -156.903779+.03))#; cachedir="/tmp")

function AOI(center::LLA, radius=1000::Number; style="osm-bright", width=1024, cachedir=nothing)
    mapimg = nothing

    radius_inDeg = radius /1852/60 # metre->degree(latitude)
    lat1 = center.lat + radius_inDeg
    lon1 = center.lon - radius_inDeg
    lat2 = center.lat - radius_inDeg
    lon2 = center.lon + radius_inDeg
    topleft = LLA(lat1, lon1); bottomright = LLA(lat2, lon2)

    @debug "horizontal: " *string(euclidean_distance(LLA(lat1, lon1), LLA(lat1, lon2)))
    @debug "vertical:   " *string(euclidean_distance(LLA(lat1, lon1), LLA(lat2, lon1)))

    height = round(Int, width * abs(lat2 - lat1) / abs(lon2 - lon1))
    try
      cachedir == nothing && (cachedir = joinpath(tempdir(), "arltoolkit", "cache"))
      filename = joinpath(cachedir, "map-$style-$lon1-$lat1-$lon2-$lat2-$width-$height.png")
      if isfile(filename)
        mapimg = FileIO.load(filename)
      else
        if "GEOAPIFY_APIKEY" ∈ keys(ENV)
          apikey = ENV["GEOAPIFY_APIKEY"]
          mkpath(cachedir)
          download("https://maps.geoapify.com/v1/staticmap?style=$style&format=png&" *
            "area=rect:$lon1,$lat1,$lon2,$lat2&width=$width&height=$height&apiKey=$apikey",
            filename)
          mapimg = FileIO.load(filename)
        end
      end
    catch ex
      @warn "Could not load map: $ex"
    end
    AOI(topleft, bottomright, mapimg)
end

AOI(LLA(20.734398, -156.903779), 1000)

function mean(lla_arr::Vector{LLA{Float64}})#AbstractVector{<:LLA})#
    pos2geo(reduce( +, geo2pos.(lla_arr, Ref(lla_arr[1]))) ./ length(lla_arr), lla_arr[1] )
end

function geo2pos(geo::LLA, origin::LLA)
  ENUfromLLA(origin, wgs84)(geo)[1:3]
  # geo.alt == 0.0 && (pos = pos[1:2])
  # pos
end
angle2pt(θ) = [sin(θ); cos(θ); 0]
angle2pt_tuple(θ) = (sin(θ), cos(θ))

## bearing zero radians from north, clockwise positive
function pts_from_location(location, bearing, len)#AbstractVector{<:LLA})#
  @debug location, bearing, len
  @debug geo2pos(location, location)
  @debug angle2pt(bearing)*len
  pos2geo( geo2pos(location, location) + angle2pt(bearing)*len , location)
end

function plot_nodes(node_loc, field_log=nothing; radius=1000)
    centre = mean(node_loc)
    aoi = AOI(centre, radius, cachedir="/tmp")
    @debug centre
    # @debug aoi
    plot(aoi, legend=nothing)
    scatter!(centre; alpha=0.1, label="centre")
    scatter!(node_loc; label="node")

    !isnothing(field_log) && annotate!(field_log.gps_longitude, field_log.gps_latitute, 
        text.(field_log.Device .*"__" .* string.(field_log.depth_m) .* "m", 8, :top, :red))
    length(node_loc) == 2 && title!("Distance: " *string(round(euclidean_distance(node_loc[1],node_loc[2]), digits=2))* "m")
    plot!()
end

function plot_arrow!(location, bearing, len, arrow_angle=π/8; kwargs...)
  @debug location, bearing, len
  plot!( [location location; pts_from_location(location, bearing + π-arrow_angle, len) pts_from_location(location, bearing  + π+arrow_angle, len)] ; kwargs...)
end

function palette_continuous(arr_scale; regular=true, colortype=:colorful)
    if regular
        scale =  range(0, 1, length(arr_scale)) #regular
    else
        scale = arr_scale .- minimum(arr_scale)
        scale = scale ./ maximum(scale)
    end
    # scale 
    
    if colortype == :colorful
      HSV.( scale*250,1,1) |> reverse
    elseif colortype == :red
      RGB.( scale,0,0)
    end
    # 
end

filter_lla(gpx, start_dts, end_dts) = filter(:time =>  <( ZonedDateTime(end_dts) ), filter(:time => >( ZonedDateTime(start_dts) ), gpx))
# seg = filter(:time =>  <( ZonedDateTime("2022-09-27T10:30:00.000-1000") ), filter(:time => >( ZonedDateTime("2022-09-27T10:28:00.000-1000") ), gpx))