from sparsemedoid import clustering
import pandas as pd
import numpy as np
from sklearn import metrics
import os

trials = []

clusters = [2, 3, 4]
distance_types = ["glower", "wishart", "podani"]
normalization_param = [1.01, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

result_labels = []

total_runs = (
    len(trials) * len(clusters) * len(distance_types) * len(normalization_param)
)

Scores = np.zeros((1, total_runs))
iter1 = 0

for trial in trials:

    # import section of data

    X = np.ones((10, 10))  # actually read in data from file

    P = X.shape[0]
    N = X.shape[1]

    runs_per_trial = total_runs / len(trials)

    all_feature_weights = np.zeros((N, runs_per_trial))
    all_cluster_labels = np.zeros((P, runs_per_trial))

    for K in clusters:

        iter2 = 0
        for distance in distance_types:

            for S in normalization_param:

                results_path_prefix = f"{trial}_K={K}_dist={distance}_S={S}"
                result_labels.append(results_path_prefix)

                print(f"N={N} | K={K} | {distance} {S} started")

                (
                    cluster_labels,
                    feature_weights,
                    feature_order,
                    weighted_distances,
                ) = clustering.sparse_kmedoids(
                    X,
                    distance_type=distance,
                    k=K,
                    s=S,
                    max_attempts=6,
                    method="pam",
                    init="build",
                    max_iter=100,
                    random_state=None,
                )

                Scores[0, iter1] += metrics.silhouette_score(
                    weighted_distances, cluster_labels, metric="precomputed"
                )

                all_feature_weights[:, iter2] = feature_weights
                all_cluster_labels[:, iter2] = cluster_labels

                iter1 += 1
                iter2 += 1

                # Save the outputs
                # weighted_distances_df = pd.DataFrame(weighted_distances)
                # weighted_distances_df.to_csv(
                #    "results/" + results_path_prefix + "_distance_matrix.csv", index=False
                # )

    feature_weights_df = pd.DataFrame(all_feature_weights)
    feature_weights_df.to_csv("results/" + trial + "_feature_weights.csv", index=False)

    cluster_labels_df = pd.DataFrame(all_cluster_labels)
    cluster_labels_df.to_csv("results/" + trial + "_cluster_labels.csv", index=False)

scores_df = pd.DataFrame(Scores)
scores_df.to_csv("results/scores.csv", index=False)
