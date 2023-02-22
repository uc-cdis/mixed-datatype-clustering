import numpy as np
import pandas as pd

from sparsemedoid import clustering


# read binary data from txt by excluding the 1st row and the 1st column
binary_data = np.loadtxt(
    "binaryData_cleaned.txt",
    dtype="str",
    delimiter="\t",
    skiprows=1,
    usecols=range(1, 572),
)

num_data = np.loadtxt(
    "propogatedData_cleaned.txt",
    dtype="float",
    delimiter="\t",
    skiprows=1,
    usecols=range(1, 514),
)

for feature in range(0, num_data.shape[1]):
    num_data[:, feature] = (num_data[:, feature] - np.min(num_data[:, feature])) / (
        np.max(num_data[:, feature] - np.min(num_data[:, feature]))
    )

X = np.concatenate((num_data, binary_data), axis=1, dtype=object)

distance = "gower"
K = 5
hyperparams = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]

cluster_labels = np.zeros((X.shape[0], len(hyperparams)))
feature_weights = np.zeros((X.shape[1], len(hyperparams)))

i = 0
# for distance in distance_types:
for S in hyperparams:

    (
        sparse_cluster_labels,
        sparse_feature_weights,
        sparse_feature_order,
    ) = clustering.sparse_kmedoids(
        X,  # original data NxP
        distance_type=distance,
        k=K,
        s=S,
        max_attempts=8,
        method="pam",
        init="build",
        max_iter=100,
        random_state=123,
    )

    cluster_labels[:, i] = sparse_cluster_labels
    feature_weights[:, i] = sparse_feature_weights
    i += 1

    print(f"Hyperparameter {S} done")

cluster_labels_df = pd.DataFrame(cluster_labels, columns=list(map(str, hyperparams)),)
cluster_labels_df.to_csv(f"sparse_cluster_labels_{distance}.csv", index=False)

feature_weights_df = pd.DataFrame(feature_weights, columns=list(map(str, hyperparams)),)
feature_weights_df.to_csv(f"sparse_feature_weights_{distance}.csv", index=False)
