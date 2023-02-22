import numpy as np
from sklearn.datasets import make_blobs
from sklearn.metrics import adjusted_rand_score
import pandas as pd
from sklearn import metrics
from sklearn.cluster import KMeans


from sparsemedoid import clustering
from sparsemedoid.distfuncs import (
    weighted_distance_matrix,
    generalized_distance_function,
)
from sparsemedoid.subfuncs import sort_datatypes, kmedoid_clusters

cluster_std = 3.0
initial_cluster_seed = 123
Scores = np.zeros((2, 45))

for cluster_trial in range(0, 50):
    print(f"Trial - {cluster_trial}")

    for N in [300]:  # [200, 300, 500]:

        for K in [3]:  # [3, 5, 7, 9]:

            cluster_seed = initial_cluster_seed + cluster_trial

            A1 = np.ones([N, 1]).astype(str)
            A2 = np.random.randint(2, size=(N, 1)).astype(str)

            B, target = make_blobs(
                n_samples=N,
                n_features=20,
                centers=K,
                cluster_std=cluster_std,
                shuffle=True,
                random_state=cluster_seed,
            )
            B1 = B[:, :10]
            B2 = B[:, 10:]
            for feat in range(0, 5):
                B2[:, feat] = 1.0 * (B2[:, feat] < np.mean(B2[:, feat]))

            B2[:, 5] = np.digitize(B2[:, 5], np.percentile(B2[:, 5], [33.3, 66.6]))
            B2[:, 6] = np.digitize(B2[:, 6], np.percentile(B2[:, 6], [25, 50, 75]))
            B2[:, 7] = np.digitize(B2[:, 7], np.percentile(B2[:, 7], [25, 50, 75]))
            B2[:, 8] = np.digitize(B2[:, 8], np.percentile(B2[:, 8], [20, 40, 60, 80]))
            B2[:, 9] = np.digitize(B2[:, 9], np.percentile(B2[:, 9], [20, 40, 60, 80]))

            B2 = B2.astype(str)

            C1 = np.random.randint(3, size=(N, 1)).astype(str)
            C2 = np.random.randint(4, size=(N, 1)).astype(str)

            X = np.concatenate((A1, A2, B1, B2, C1, C2), axis=1, dtype=object)

            hyperparams = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
            distance_types = ["gower", "wishart", "podani", "huang"]

            trials = len(distance_types) * (len(hyperparams) + 1) + 1

            # Scores = np.zeros((2, trials))
            all_cluster_labels = np.zeros((X.shape[0], trials))
            all_feature_weights = np.zeros((X.shape[1], trials))

            iter = 0
            for distance in distance_types:

                for S in hyperparams:
                    print(f"N={N} | K={K} | {distance} {S} started")

                    (
                        cluster_labels,
                        feature_weights,
                        feature_order1,
                        weighted_distances,
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

                    Scores[0, iter] += metrics.silhouette_score(
                        weighted_distances, cluster_labels, metric="precomputed"
                    )
                    Scores[1, iter] += adjusted_rand_score(target, cluster_labels)

                    # all_feature_weights[:, iter] = feature_weights
                    # all_cluster_labels[:, iter] = cluster_labels

                    weighted_distances_df = pd.DataFrame(weighted_distances)
                    # weighted_distances_df.to_csv(f"generated_data/sparse/cluster_std_{cluster_std}/{N}/{K}/{distance}_{S}_distance_matrix.csv", index=False)
                    iter += 1

                for unweighted in ["unweighted"]:
                    print(f"N={N} | K={K} | {distance} {unweighted} started")

                    n, p = X.shape

                    x_numeric, x_binary, x_categoric, feature_order2 = sort_datatypes(
                        X, p
                    )

                    feature_counts = {
                        "Numeric": x_numeric.shape[1],
                        "Binary": x_binary.shape[1],
                        "Categoric": x_categoric.shape[1],
                    }

                    per_feature_distances = generalized_distance_function(
                        x_numeric, x_binary, x_categoric, distance, feature_counts
                    )

                    unweighted_feature_weights = np.ones(p) * (1 / np.sqrt(p))
                    unweighted_weighted_distances = weighted_distance_matrix(
                        per_feature_distances,
                        unweighted_feature_weights,
                        distance,
                        feature_counts,
                    )

                    unweighted_weighted_distances_df = pd.DataFrame(
                        unweighted_weighted_distances
                    )
                    # unweighted_weighted_distances_df.to_csv(
                    #    f"generated_data/sparse/cluster_std_{cluster_std}/{N}/{K}/unweighted_{distance}_distance_matrix.csv", index=False
                    # )

                    norm_cluster_labels = kmedoid_clusters(
                        unweighted_weighted_distances,
                        k=K,
                        method="pam",
                        init="build",
                        max_iter=100,
                        random_state=None,
                    )

                    Scores[0, iter] += metrics.silhouette_score(
                        unweighted_weighted_distances,
                        norm_cluster_labels,
                        metric="precomputed",
                    )

                    Scores[1, iter] += adjusted_rand_score(target, norm_cluster_labels)

                    # all_feature_weights[:, iter] = unweighted_feature_weights
                    # all_cluster_labels[:, iter] = norm_cluster_labels

                    iter += 1

            for k_means in ["kmeans"]:
                print(f"N={N} | K={K} | {k_means} started")
                col_ints = list(range(1, 25))
                col_names = [str(col) for col in col_ints]
                df_X = pd.DataFrame(X, columns=col_names)
                encode = [
                    "1",
                    "2",
                    "13",
                    "14",
                    "15",
                    "16",
                    "17",
                    "18",
                    "19",
                    "20",
                    "21",
                    "22",
                    "23",
                    "24",
                ]
                one_hot_encodings = pd.get_dummies(df_X[encode])
                df_X = df_X.drop(encode, axis=1)
                df_X = df_X.join(one_hot_encodings)
                encoded_X = df_X.to_numpy()

                kmeans_clusters = KMeans(n_clusters=K).fit(encoded_X)
                kmeans_labels = kmeans_clusters.labels_

                Scores[0, iter] += metrics.silhouette_score(encoded_X, kmeans_labels)
                Scores[1, iter] += adjusted_rand_score(target, kmeans_labels)

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
                "huang_unweighted",
                "kmeans_unweighted",
            ]

            labels_df = pd.DataFrame(all_cluster_labels, columns=cols)
            weights_df = pd.DataFrame(all_feature_weights, columns=cols)

Scores = Scores / 50
Scores_df = pd.DataFrame(Scores, columns=cols)

Scores_df.to_csv(
    f"generated_data/sparse/cluster_std_{cluster_std}/{N}/{K}/mean_generated_data_scores.csv",
    index=False,
)
# Scores_df.to_csv(f"generated_data/sparse/cluster_std_{cluster_std}/{N}/{K}/generated_data_scores.csv", index=False)
# labels_df.to_csv(f"generated_data/sparse/cluster_std_{cluster_std}/{N}/{K}/generated_data_labels.csv", index=False)
# weights_df.to_csv(f"generated_data/sparse/cluster_std_{cluster_std}/{N}/{K}/generated_data_weights.csv", index=False)
