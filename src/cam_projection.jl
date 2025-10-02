using ImageProjectiveGeometry
using JLD2
using LinearAlgebra
using JSON

export get_cam_fov_fromparam

function Camera2(camera_params)
    Camera(
        fx  = camera_params["mtx"][1,1],
        fy  = camera_params["mtx"][2,2],
        ppx = camera_params["mtx"][3,1],
        ppy = camera_params["mtx"][3,2],
        k1  = camera_params["dist"][1],
        k2  = camera_params["dist"][2],
        k3  = camera_params["dist"][5],
        p1  = camera_params["dist"][3],
        p2  = camera_params["dist"][4],
        skew= camera_params["mtx"][2,1],
        cols= camera_params["imsize"][1],
        rows= camera_params["imsize"][2]
    )
end

vect_angle(a,b) = acos(clamp(a⋅b/(norm(a)*norm(b)), -1, 1))
vect_angled(a,b) = acosd(clamp(a⋅b/(norm(a)*norm(b)), -1, 1))

# sph2cart(r, θ, ϕ) = r. [; ; ] 

function get_cam_fov_fromparam(fname)
    if endswith(fname, ".json")
        camera_params_file = JSON.parsefile(fname)
        camera_params = Dict(
            "mtx" => hcat(camera_params_file["pinhole"]["camera_matrix"]...),
            "dist" => camera_params_file["pinhole"]["dist_coeffs"][1],
            "imsize" => (camera_params_file["resolution"]["width"], camera_params_file["resolution"]["height"])
        )
    elseif endswith(fname, ".h5py") || endswith(fname, ".h5")
        camera_params = load(fname) #"/Users/abel/Documents/data_res/aspod/cam_calib/aspod2/Vid_20131219_105014_s30.h5")
    else
        error("Unsupported file format: $fname")
    end
    
    # return camera_params
    
    cam = Camera2(camera_params)

    horizontal_angle = vect_angled( imagept2ray(cam, cam.rows/2, 1), imagept2ray(cam, cam.rows/2, cam.cols) ) # 62.61721188568244
    vertical_angle = vect_angled( imagept2ray(cam, 1, cam.cols/2), imagept2ray(cam, cam.rows, cam.cols/2) ) # 35.793211268714096
    diagonal_angle = vect_angled( imagept2ray(cam, 1, 1), imagept2ray(cam, cam.rows, cam.cols)) # 71.6855447884958

    return (;horizontal_angle, vertical_angle, diagonal_angle)    
end

"""
    get_cam_fov_fromMeasuredDistancesofFOV(d_horizontal, d_vertical, d_diagonal=missing, obj_distance=1.0)

Calculate the camera's field of view (FOV) angles (horizontal, vertical, and diagonal) in radians, given the measured distances of the FOV and the distance to the object.

# Arguments
- `d_horizontal::Number`: Measured horizontal distance of the FOV.
- `d_vertical::Number`: Measured vertical distance of the FOV.
- `d_diagonal::Number=missing`: Measured diagonal distance of the FOV (optional).
- `obj_distance::Number=1.0`: Distance from the camera to the object (default is 1.0).

# Returns
Named tuple with:
- `horizontal_angle`: Horizontal FOV angle in radians.
- `vertical_angle`: Vertical FOV angle in radians.
- `diagonal_angle`: Diagonal FOV angle in radians (or `missing` if not provided).
"""

"""
    get_cam_fov_fromMeasuredDistancesofFOVd2(args...)

Calculate the camera's field of view (FOV) angles in degrees, given the measured distances of the FOV and the distance to the object.

# Arguments
- `args...`: Arguments passed to `get_cam_fov_fromMeasuredDistancesofFOV`.

# Returns
Named tuple with:
- `horizontal_angle`: Horizontal FOV angle in degrees.
- `vertical_angle`: Vertical FOV angle in degrees.
- `diagonal_angle`: Diagonal FOV angle in degrees (or `missing` if not provided).
"""
get_cam_fov_fromMeasuredDistancesofFOV(d_horizontal, d_vertical, d_diagonal=missing, obj_distance=1.0) = (
    horizontal_angle = 2.0 * atan(d_horizontal/(2.0*obj_distance)),
    vertical_angle = 2.0 * atan(d_vertical/(2.0*obj_distance)),
    diagonal_angle = 2.0 * atan(d_diagonal/(2.0*obj_distance)),
)


"""
    get_cam_fov_fromMeasuredDistancesofFOVd2(args...)

Calculate the camera's field of view (FOV) angles in degrees, given the measured distances of the FOV and the distance to the object.

# Arguments
- `args...`: Arguments passed to `get_cam_fov_fromMeasuredDistancesofFOV`.

# Returns
Named tuple with:
- `horizontal_angle`: Horizontal FOV angle in degrees.
- `vertical_angle`: Vertical FOV angle in degrees.
- `diagonal_angle`: Diagonal FOV angle in degrees (or `missing` if not provided).

"""
get_cam_fov_fromMeasuredDistancesofFOVd(args...) = map(x-> x isa Number ? rad2deg(x) : x, get_cam_fov_fromMeasuredDistancesofFOV(args...))
  

    

