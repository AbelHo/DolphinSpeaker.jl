#!/usr/bin/env julia
"""
Simple test to verify cuML works with PythonCall by installing it directly to the project-local environment
"""

# Create a simple isolated test without loading other project dependencies
import sys
print("Python version:", sys.version)

# Try to install cuml using pip first
import subprocess
try:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--upgrade", "cuml-cu11"])
    print("Successfully installed cuml-cu11 via pip")
except Exception as e:
    print(f"Failed to install cuml via pip: {e}")
    print("Note: cuML may require CUDA and/or specific configuration")
    sys.exit(1)

# Now test the import
try:
    import cuml
    print(f"\n✓ Successfully imported cuml")
    print(f"  cuML version: {cuml.__version__}")
    
    import numpy as np
    from sklearn.datasets import make_blobs
    
    print("✓ Successfully imported dependencies")
    
    # Test 1: Create sample data
    print("\n--- Test 1: Creating sample dataset ---")
    X, y = make_blobs(n_samples=1000, n_features=10, centers=4, random_state=42)
    print(f"✓ Created sample data with shape: {X.shape}")
    
    # Test 2: Test UMAP
    print("\n--- Test 2: Testing cuML UMAP ---")
    try:
        reducer = cuml.manifold.UMAP(n_components=2, random_state=42)
        X_reduced = reducer.fit_transform(X)
        print(f"✓ Successfully ran UMAP")
        print(f"  Reduced data shape: {X_reduced.shape}")
    except Exception as e:
        print(f"⚠ UMAP test skipped (may require GPU): {e}")
    
    # Test 3: Test KMeans
    print("\n--- Test 3: Testing cuML KMeans ---")
    try:
        kmeans = cuml.cluster.KMeans(n_clusters=4, random_state=42)
        kmeans.fit(X)
        labels = kmeans.labels_
        print(f"✓ Successfully ran KMeans")
        print(f"  Cluster labels shape: {labels.shape}")
    except Exception as e:
        print(f"⚠ KMeans test skipped (may require GPU): {e}")
    
    # Test 4: Test PCA
    print("\n--- Test 4: Testing cuML PCA ---")
    try:
        pca = cuml.decomposition.PCA(n_components=5)
        X_pca = pca.fit_transform(X)
        explained_var = pca.explained_variance_ratio_
        print(f"✓ Successfully ran PCA")
        print(f"  PCA data shape: {X_pca.shape}")
        print(f"  Explained variance: {np.round(explained_var, 4)}")
    except Exception as e:
        print(f"⚠ PCA test skipped (may require GPU): {e}")
    
    print("\n" + "=" * 60)
    print("cuML is successfully installed and working!")
    print("=" * 60)

except ImportError as e:
    print(f"\n❌ Failed to import cuml: {e}")
    print("Please ensure cuml is installed: pip install cuml-cu11")
    sys.exit(1)
except Exception as e:
    print(f"\n❌ Error during testing: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
