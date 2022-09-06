import numpy as np
from sklearn.datasets import make_blobs
from sklearn.metrics import adjusted_rand_score
import pandas as pd
from sklearn import metrics


from sparsemedoid import clustering
from sparsemedoid import subfuncs, distfuncs

n = 200
K = 5


B, target = make_blobs(
    n_samples=n,
    n_features=20,
    centers=K,
    cluster_std=1.7,
    shuffle=True,
    random_state=123,
)
B1 = B[:, :10]
B2 = B[:, 10:]
for feat in range(0, 10):

    B2[:, feat] = 1.0 * (B2[:, feat] < np.mean(B2[:, feat]))

B2 = B2.astype(str)

A1 = np.ones([n, 1]).astype(str)

A2 = np.random.randint(2, size=(n, 1)).astype(str)

C1 = np.random.randint(3, size=(n, 1)).astype(str)

C2 = np.random.randint(4, size=(n, 1)).astype(str)

X = np.concatenate((A1, A2, B1, B2, C1, C2), axis=1, dtype=object)


Scores = np.zeros((2, 48))
all_cluster_labels = np.zeros((X.shape[0], 48))
all_feature_weights = np.zeros((X.shape[1], 48))

distance_types = ["gower", "wishart", "podani", "huang"]
iter = 0
for distance in distance_types:
    print(f"{distance} started")

    hyperparams = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]

    for S in hyperparams:

        (
            sparse_cluster_labels,
            sparse_feature_weights,
            weight_difference1,
            feature_order1,
            sparse_weighted_distances,
        ) = clustering.sparse_kmedoids(
            X,  # original data NxP
            distance_type=distance,
            k=K,
            s=S,
            max_attempts=6,
            method="pam",
            init="build",
            max_iter=100,
            random_state=None,
        )

        Scores[0, iter] = metrics.silhouette_score(
            sparse_weighted_distances, sparse_cluster_labels, metric="precomputed"
        )

        Scores[1, iter] = adjusted_rand_score(target, sparse_cluster_labels)
        all_feature_weights[:, iter] = sparse_feature_weights
        all_cluster_labels[:, iter] = sparse_cluster_labels

        weighted_distances_df = pd.DataFrame(sparse_weighted_distances)
        weighted_distances_df.to_csv(
            f"generated_data/sparse_{distance}_distance_matrix_{S}.csv", index=False
        )

        iter += 1
        print(f"\t {distance} {S} done")

    (
        spectral_cluster_labels,
        spectral_feature_weights,
        feature_order2,
        spectral_weighted_distances,
    ) = clustering.spectral_kmedoids(
        X,
        distance_type=distance,
        k=K,
        method="pam",
        init="build",
        max_iter=100,
        random_state=None,
    )

    Scores[0, iter] = metrics.silhouette_score(
        spectral_weighted_distances, spectral_cluster_labels, metric="precomputed"
    )

    Scores[1, iter] = adjusted_rand_score(target, spectral_cluster_labels)

    all_feature_weights[:, iter] = spectral_feature_weights
    all_cluster_labels[:, iter] = spectral_cluster_labels

    weighted_distances_df = pd.DataFrame(spectral_weighted_distances)
    weighted_distances_df.to_csv(
        f"generated_data/spectral_{distance}_distance_matrix.csv", index=False
    )
    print(f"\t {distance} spectral done")

    n, p = X.shape
    x_numeric, x_binary, x_categoric, feature_order3 = subfuncs.sort_datatypes(X, p)
    feature_counts = {
        "Numeric": x_numeric.shape[1],
        "Binary": x_binary.shape[1],
        "Categoric": x_categoric.shape[1],
    }

    per_feature_distances3 = distfuncs.generalized_distance_function(
        x_numeric, x_binary, x_categoric, distance, feature_counts
    )

    norm_feature_weights = np.ones(p) * (1 / np.sqrt(p))
    norm_weighted_distances = distfuncs.weighted_distance_matrix(
        per_feature_distances3, norm_feature_weights, distance, feature_counts
    )

    weighted_distances_df = pd.DataFrame(norm_weighted_distances)
    weighted_distances_df.to_csv(
        f"generated_data/unweighted_{distance}_distance_matrix.csv", index=False
    )

    norm_cluster_labels = subfuncs.kmedoid_clusters(
        norm_weighted_distances,
        k=5,
        method="pam",
        init="build",
        max_iter=100,
        random_state=None,
    )

    Scores[0, iter] = metrics.silhouette_score(
        norm_weighted_distances, norm_cluster_labels, metric="precomputed"
    )

    Scores[1, iter] = adjusted_rand_score(target, norm_cluster_labels)

    all_feature_weights[:, iter] = norm_feature_weights
    all_cluster_labels[:, iter] = norm_cluster_labels

    iter += 1
    print(f"\t {distance} unweighted done")


cols = [
    "gower_sparse_1.2",
    "gower_sparse_1.4",
    "gower_sparse_1.6",
    "gower_sparse_1.8",
    "gower_sparse_2.0",
    "gower_sparse_2.2",
    "gower_sparse_2.4",
    "gower_sparse_2.6",
    "gower_sparse_2.8",
    "gower_sparse_3.0",
    "gower_spectral",
    "gower_unweighted",
    "wishart_sparse_1.2",
    "wishart_sparse_1.4",
    "wishart_sparse_1.6",
    "wishart_sparse_1.8",
    "wishart_sparse_2.0",
    "wishart_sparse_2.2",
    "wishart_sparse_2.4",
    "wishart_sparse_2.6",
    "wishart_sparse_2.8",
    "wishart_sparse_3.0",
    "wishart_spectral",
    "wishart_unweighted",
    "podani_sparse_1.2",
    "podani_sparse_1.4",
    "podani_sparse_1.6",
    "podani_sparse_1.8",
    "podani_sparse_2.0",
    "podani_sparse_2.2",
    "podani_sparse_2.4",
    "podani_sparse_2.6",
    "podani_sparse_2.8",
    "podani_sparse_3.0",
    "podani_spectral",
    "podani_unweighted",
    "huang_sparse_1.2",
    "huang_sparse_1.4",
    "huang_sparse_1.6",
    "huang_sparse_1.8",
    "huang_sparse_2.0",
    "huang_sparse_2.2",
    "huang_sparse_2.4",
    "huang_sparse_2.6",
    "huang_sparse_2.8",
    "huang_sparse_3.0",
    "huang_spectral",
    "huang_unweighted",
]

Scores_df = pd.DataFrame(Scores, columns=cols)
labels_df = pd.DataFrame(all_cluster_labels, columns=cols)
weights_df = pd.DataFrame(all_feature_weights, columns=cols)

Scores_df.to_csv(
    "generated_data/generated_data_silhouette_and_ARI_scores.csv", index=False
)
labels_df.to_csv("generated_data/generated_data_labels.csv", index=False)
weights_df.to_csv("generated_data/generated_data_weights.csv", index=False)
