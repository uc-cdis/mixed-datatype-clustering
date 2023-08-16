import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from sklearn_extra.cluster import KMedoids
from sklearn.datasets import make_blobs
import pandas as pd
from sklearn import metrics
from sklearn.cluster import KMeans
import os
import sys

# get absolute path to the directory where the script resides
current_dir = os.path.dirname(os.path.abspath(__file__))

# get absolute path to the parent directory
parent_dir = os.path.dirname(current_dir)

# add the parent directory to the sys.path
sys.path.append(parent_dir)

from sparsemedoid import clustering

# Keep it as a list for now so if we want to change back to iterating
# within the script it is easy to do so
hyperparams = [float(sys.argv[1])]
distance_types = [sys.argv[2]]
clusters = [int(sys.argv[3])]
maf_file = sys.argv[4]

maf_df = pd.read_csv(maf_file, index_col=False)

X_full = maf_df.values
X = X_full

N = X.shape[0]
P = X.shape[1]

# trials = len(clusters) * len(distance_types) * len(hyperparams)

Scores = np.zeros((1, 1))
all_feature_weights = np.zeros((P, 1))
all_feature_orders = np.zeros((P, 1))
all_cluster_labels = np.zeros((N, 1))

result_labels = []

for K in clusters:

    iter = 0
    for distance in distance_types:

        for S in hyperparams:

            results_path_prefix = f"P={P}_K={K}_dist={distance}_S={S}"
            result_labels.append(results_path_prefix)

            print(f"P={P} | K={K} | {distance} {S} started")

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

            Scores[0, iter] += metrics.silhouette_score(
                weighted_distances, cluster_labels, metric="precomputed"
            )

            all_feature_weights[:, iter] = feature_weights
            all_cluster_labels[:, iter] = cluster_labels
            feature_order_vec = (
                feature_order["Numerical Features"]
                + feature_order["Binary Features"]
                + feature_order["Categorical Features"]
            )
            all_feature_orders[:, iter] = feature_order_vec

            weighted_distances_df = pd.DataFrame(weighted_distances)
            os.makedirs(
                os.path.dirname(
                    "results/" + results_path_prefix + "_distance_matrix.csv"
                ),
                exist_ok=True,
            )
            weighted_distances_df.to_csv(
                "results/" + results_path_prefix + "_distance_matrix.csv", index=False
            )

            iter += 1

labels_df = pd.DataFrame(all_cluster_labels, columns=result_labels)
os.makedirs(
    os.path.dirname("results/" + results_path_prefix + "_cluster_labels.csv"),
    exist_ok=True,
)
labels_df.to_csv("results/" + results_path_prefix + "_cluster_labels.csv", index=False)

weights_df = pd.DataFrame(all_feature_weights, columns=result_labels)
os.makedirs(
    os.path.dirname("results/" + results_path_prefix + "_cluster_weights.csv"),
    exist_ok=True,
)
weights_df.to_csv(
    "results/" + results_path_prefix + "_cluster_weights.csv", index=False
)

Scores_df = pd.DataFrame(Scores, columns=result_labels)
os.makedirs(
    os.path.dirname("results/" + results_path_prefix + "_cluster_scores.csv"),
    exist_ok=True,
)
Scores_df.to_csv("results/" + results_path_prefix + "_cluster_scores.csv", index=False)

orders_df = pd.DataFrame(all_feature_orders, columns=result_labels)
os.makedirs(
    os.path.dirname("results/" + results_path_prefix + "_feature_orders.csv"),
    exist_ok=True,
)
orders_df.to_csv("results/" + results_path_prefix + "_feature_orders.csv", index=False)
