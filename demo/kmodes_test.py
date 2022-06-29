import numpy as np

from sklearn.datasets import make_blobs
from sklearn.metrics import adjusted_rand_score

from kmodes.kprototypes import KPrototypes


n = 50
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


kprototype = KPrototypes(n_jobs=1, n_clusters=K, init="Huang", random_state=0)

kprototype.fit(X=X, categorical=[0, 5], sample_weight=[1, 1, 1, 1, 1, 1])

labels = kprototype.fit_predict(X, categorical=[0, 5])


print("Sparse ARI = ", adjusted_rand_score(target, labels))
