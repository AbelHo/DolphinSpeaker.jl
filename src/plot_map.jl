import PlotlyJS
function maps1()
    marker = PlotlyJS.attr(size=[20, 30, 15, 10],
                    color=[10, 20, 40, 50],
                    cmin=0,
                    cmax=50,
                    colorscale="Greens",
                    colorbar=attr(title="Some rate",
                                ticksuffix="%",
                                showticksuffix="last"),
                    line_color="black")
    trace = PlotlyJS.scattergeo(;mode="markers+lines", locations=["FRA", "DEU", "RUS", "ESP"],
                        marker=marker, name="Europe Data")
    layout = PlotlyJS.Layout(geo_scope="europe", geo_resolution=50, width=500, height=550,
                    margin=attr(l=0, r=0, t=10, b=0))
    PlotlyJS.plot(trace, layout)
end

function maps2()
    marker = PlotlyJS.attr(size=[20, 30, 15, 10],
                    color=[10, 20, 40, 50],
                    cmin=0,
                    cmax=50,
                    colorscale="Greens",
                    colorbar=attr(title="Some rate",
                                ticksuffix="%",
                                showticksuffix="last"),
                    line_color="black")
    trace = PlotlyJS.scattergeo(;mode="markers+lines", locations=[[20.734398, -156.903779], [20.734398, -156.903779],[20.734398, -156.903779],[20.734398, -156.903779]],
                        marker=marker, name="Europe Data")
    layout = PlotlyJS.Layout(geo_center_lat=20.734398, geo_center_lon=-156.903779, geo_resolution=50, width=500, height=550,
                    margin=attr(l=0, r=0, t=10, b=0))
    PlotlyJS.plot(trace, layout)
end

using GeoMakie, CairoMakie
using GeoMakie.GeoJSON
function map_makie()
    # First, make a surface plot
    lons = -160:-140
    lats = 10:30
    field = [exp(cosd(l)) + 3(y/90) for l in lons, y in lats]

    fig = GeoMakie.Figure()
    ax = GeoMakie.GeoAxis(fig[1,1])
    sf = GeoMakie.surface!(ax, lons, lats, field; shading = false)
    cb1 = GeoMakie.Colorbar(fig[1,2], sf; label = "field", height = Relative(0.65))

    using GeoMakie.GeoJSON
    countries_file = download("https://datahub.io/core/geo-countries/r/countries.geojson")
    countries = GeoJSON.read(read(countries_file, String))

    n = length(countries)
    hm = poly!(ax, countries; color= 1:n, colormap = :dense,
        strokecolor = :black, strokewidth = 0.5,
    )
    translate!(hm, 0, 0, 100) # move above surface plot

    fig
end