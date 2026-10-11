#!/usr/bin/env julia
"""
Test script to verify cuml installation and functionality with PythonCall.jl
"""

using PythonCall

divider = repeat("=", 60)
println(divider)
println("Testing cuML with PythonCall.jl")
println(divider)

try
    # Import cuml
    cuml = pyimport("cuml")
    println("\n✓ Successfully imported cuml")
    println("  cuML version: $(pyconvert(String, cuml.__version__))")
    
    # Import additional required packages
    np = pyimport("numpy")
    sklearn = pyimport("sklearn.datasets")
    
    println("✓ Successfully imported numpy and sklearn")
    
    # Test 1: Create sample data using sklearn
    println("\n--- Test 1: Creating sample dataset ---")
    X, y = sklearn.make_blobs(
        n_samples=1000,
        n_features=10,
        centers=4,
        random_state=42
    )
    shape = pyconvert(Tuple, X.shape)
    println("✓ Created sample data with shape: $shape")
    
    # Test 2: Test UMAP (a common cuML algorithm)
    println("\n--- Test 2: Testing cuML UMAP ---")
    reducer = cuml.manifold.UMAP(
        n_components=2,
        random_state=42
    )
    X_reduced = reducer.fit_transform(X)
    reduced_shape = pyconvert(Tuple, X_reduced.shape)
    println("✓ Successfully ran UMAP")
    println("  Reduced data shape: $reduced_shape")
    
    # Test 3: Test KMeans clustering
    println("\n--- Test 3: Testing cuML KMeans ---")
    kmeans = cuml.cluster.KMeans(
        n_clusters=4,
        random_state=42
    )
    kmeans.fit(X)
    labels = kmeans.labels_
    labels_shape = pyconvert(Tuple, labels.shape)
    println("✓ Successfully ran KMeans")
    println("  Cluster labels shape: $labels_shape")
    
    # Test 4: Test PCA
    println("\n--- Test 4: Testing cuML PCA ---")
    pca = cuml.decomposition.PCA(
        n_components=5
    )
    X_pca = pca.fit_transform(X)
    pca_shape = pyconvert(Tuple, X_pca.shape)
    explained_var = pyconvert(Vector, pca.explained_variance_ratio_)
    println("✓ Successfully ran PCA")
    println("  PCA data shape: $pca_shape")
    println("  Explained variance: $(round.(explained_var, digits=4))")
    
    println("\n" * repeat("=", 60))
    println("All tests passed! ✓")
    println("cuML is successfully installed and working with PythonCall.jl")
    println(repeat("=", 60))
    
catch e
    println("\n❌ Error during testing:")
    println("  $(sprint(showerror, e))")
    rethrow()
end
