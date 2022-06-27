import numpy as np
from sklearn.datasets import make_blobs
from sklearn.metrics import adjusted_rand_score

from sparsemedoid import clustering

n = 500
K = 5

A = np.ones([n, 1]).astype(str)

B, target = make_blobs(
    n_samples=n, n_features=4, centers=K, cluster_std=0.2, shuffle=True
)


C = np.ones([n, 1]).astype(str)

for i in range(0, n):
    if i % 3 == 0:
        C[i, 0] = "0"
    elif i % 3 == 1:
        C[i, 0] = "1"
    else:
        C[i, 0] = "2"

X = np.concatenate((A, B, C), axis=1, dtype=object)


(
    cluster_labels1,
    weights1,
    weight_difference1,
    feature_order1,
) = clustering.sparse_kmedoids(
    X,
    distance_type="gower",
    k=K,
    s=1.2,
    max_attempts=6,
    method="pam",
    init="build",
    max_iter=100,
    random_state=None,
)


cluster_labels2, weights2, feature_order2 = clustering.spectral_kmedoids(
    X,
    distance_type="gower",
    k=K,
    method="pam",
    init="build",
    max_iter=100,
    random_state=None,
)


print("Sparse ARI = ", adjusted_rand_score(target, cluster_labels1))
print("Spectral ARI = ", adjusted_rand_score(target, cluster_labels2))
print("Sparse Feature weights = ", weights1)
print("Spectral Feature weights = ", weights2)
print("Sparse Feature Order = ", feature_order1)
print("Spectral Feature Order = ", feature_order2)

# print('Entropy term = ', weight_difference)
