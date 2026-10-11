# cuML Installation Guide with CondaPkg.jl and PythonCall.jl

## Summary

cuML has been added to `CondaPkg.toml` and is configured to work with PythonCall.jl. However, cuML installation requires careful management of Conda channels due to dependency conflicts.

## Current Configuration

Your `CondaPkg.toml` now includes:

```toml
channels = ["conda-forge"]

[deps]
pandas = ""
matplotlib = ""
hdbscan = ""
cuml = ""
```

##Installation Approaches

### Option 1: Using Conda Directly (Recommended if you have GPU)

Install cuML outside the Julia environment using conda/mamba with RAPIDS channels:

```bash
# Using conda with RAPIDS channels
conda create -n rapids_env -c rapidsai -c nvidia -c conda-forge cuml cuda-version=12.0

# Or using mamba (faster)
mamba create -n rapids_env -c rapidsai -c nvidia -c conda-forge cuml

# Then activate and use Julia
conda activate rapids_env
julia --project=.
```

Then in Julia, point PythonCall to use this conda environment:

```julia
using PythonCall
PythonCall.set_python_executable("/path/to/rapids_env/bin/python")
```

### Option 2: CPU-only Installation (if no GPU available)

Install the CPU-only version via pip:

```bash
pip install cuml-cpu
```

Then in Julia, cuML will be available through PythonCall.

### Option 3: Use CondaPkg with Workarounds

The challenge with CondaPkg is that MLJScikitLearnInterface requires a specific scikit-learn version from conda-forge, which conflicts with RAPIDS packages. To work around this:

1. **Temporarily remove MLJScikitLearnInterface** from your Project.toml if you don't need it
2. **Use a dedicated CondaPkg environment** by modifying Julia's environment selection

## Testing cuML Installation

Once installed, test with this Julia code:

```julia
using PythonCall

# Import cuML
cuml = pyimport("cuml")
println("cuML version: $(cuml.__version__)")

# Simple test: UMAP dimensionality reduction
import numpy as np
from sklearn.datasets import make_blobs

X, y = make_blobs(n_samples=1000, n_features=10, centers=4, random_state=42)
reducer = cuml.manifold.UMAP(n_components=2)
X_reduced = reducer.fit_transform(X)
println("UMAP output shape: $(X_reduced.shape)")
```

## Troubleshooting

###  Dependency Conflicts
If you see "strict channel priority" errors, the issue is that conda-forge packages have different version constraints than RAPIDS. Solutions:

- Use a separate conda environment dedicated to RAPIDS (Option 1)
- Use the CPU-only version (Option 2)
- Remove conflicting packages from your Project.toml temporarily

### GPU CUDA Compatibility
cuML requires CUDA 11.0 or higher. Check your CUDA version:

```bash
nvidia-smi
nvcc --version
```

###  Module Not Found Errors
If PythonCall can't find cuML, ensure:

1. cuML is installed in the Python environment PythonCall is using
2. Check which Python PythonCall is using:
   ```julia
   using PythonCall
   PythonCall.python_executable()
   ```

## Files Created

- `test_cuml.jl` - Julia test script with syntax-corrected string operations
- `test_cuml_direct.py` - Standalone Python test script
- `CUML_SETUP.md` - This documentation

## Next Steps

1. Choose an installation approach above
2. Install cuML using one of the methods
3. Run the test script to verify it works
4. Import and use cuML in your Julia code through PythonCall

For more information, see:
- [cuML Documentation](https://docs.rapids.ai/api/cuml/stable/)
- [RAPIDS Installation Guide](https://rapids.ai/install)
- [PythonCall.jl Documentation](https://github.com/cjdoris/PythonCall.jl)
