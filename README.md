# Dolphin Speaker Localization
<!--- 
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaci.github.io/PkgTemplates.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaci.github.io/PkgTemplates.jl/dev)
[![CI](https://github.com/JuliaCI/PkgTemplates.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/JuliaCI/PkgTemplates.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Codecov](https://codecov.io/gh/JuliaCI/PkgTemplates.jl/branch/master/graph/badge.svg?token=WsGRSymBmZ)](https://codecov.io/gh/JuliaCI/PkgTemplates.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

**PkgTemplates creates new Julia packages in an easy, repeatable, and customizable way.**
--->
## Installation

Install with the Julia package manager [Pkg](https://pkgdocs.julialang.org/), just like any other registered Julia package:
1. Press "]"
2.
```
add https://github.com/AbelHo/DolphinSpeaker.jl
```

## Usage
Convert .mat to .flac audio file
Format:
```
mat2flac(input_folderpath; outfilepath=output_folderpath)
```
Example:
```
mat2flac("/Users/abel/Documents/data/calf/Clicktest/20190913"; outfilepath="/Users/abel/Documents/data/calf/Clicktest/20190913/flac")
```

Process and overlay one set of video and acoustic recording, provided to results directory/folder
format:
```
process_one_set(vidfname, aufname, res_dir)
```
Example:
```
process_one_set("/Users/abel/Documents/data/calf/Clicktest/20231129/20231129_15.20.54_log.mkv", "/Users/abel/Documents/data/calf/Clicktest/20231129/20231129_15.20.54_log.flac", "/Users/abel/Documents/data_res/calf/temp/test_laptop")
```

### 

## Contributing


