

function angle2px(ang, fov_angled=[62.61721188568244,35.793211268714096,71.6855447884958], imsize=(2160, 3840))
	angle2pxd( rad2deg.(ang), fov_angled, imsize)
end

function angle2pxd(angd, fov_angled=[62.61721188568244,35.793211268714096,71.6855447884958], imsize=(2160, 3840))
	angd ./ [fov_angled[1]/2 -fov_angled[2]/2] .* [imsize[2]/2 imsize[1]/2] .+ [imsize[2]/2 imsize[1]/2]
end

p_pixels = angle2px(ang, fov_angle)
pind_vidframes = round.(Int, pind_good_inS * get_fps(vidfname)) .+ 1