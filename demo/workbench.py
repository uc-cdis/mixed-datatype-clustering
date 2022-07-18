import numpy as np
import pandas as pd
from sklearn import metrics

from sparsemedoid import subfuncs, distfuncs


# read binary data from txt by excluding the 1st row and the 1st column
binary_data = np.loadtxt(
    "binaryData_cleaned.txt",
    dtype="str",
    delimiter="\t",
    skiprows=1,
    usecols=range(1, 572),
)
# print(binary_data.shape)


num_data = np.loadtxt(
    "propogatedData_cleaned.txt",
    dtype="float",
    delimiter="\t",
    skiprows=1,
    usecols=range(1, 514),
)
# print(num_data.shape)
# print(num_data)

for feature in range(0, num_data.shape[1]):
    num_data[:, feature] = (num_data[:, feature] - np.min(num_data[:, feature])) / (
        np.max(num_data[:, feature] - np.min(num_data[:, feature]))
    )

X = np.concatenate((num_data, binary_data), axis=1, dtype=object)
# distance_types = ["gower", "wishart", "podani", "huang"]
distance = "wishart"
K = 5

hyperparams = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]

feature_weights_df = pd.read_csv(f"sparse_feature_weights_{distance}.csv")
labels_df = pd.read_csv(f"sparse_cluster_labels_{distance}.csv")

all_feature_weights = feature_weights_df.to_numpy()
all_labels = labels_df.to_numpy()

Silh_scores = np.zeros(len(hyperparams))

n, p = X.shape

x_numeric, x_binary, x_categoric, feature_order = subfuncs.sort_datatypes(X, p)

feature_counts = {
    "Numeric": x_numeric.shape[1],
    "Binary": x_binary.shape[1],
    "Categoric": x_categoric.shape[1],
}
print("NxNxP started")
per_feature_distances = distfuncs.generalized_distance_function(
    x_numeric, x_binary, x_categoric, distance, feature_counts
)
print("NxNxP done")

i = 0
for S in hyperparams:

    feature_weights = all_feature_weights[:, i]
    labels = all_labels[:, i]

    weighted_distances = distfuncs.weighted_distance_matrix(
        per_feature_distances, feature_weights, distance, feature_counts
    )

    weighted_distances_df = pd.DataFrame(weighted_distances)
    weighted_distances_df.to_csv(
        f"sparse_cluster_weighted_{distance}_{S}.csv", index=False
    )

    Silh_scores[i] = metrics.silhouette_score(
        weighted_distances, labels, metric="precomputed"
    )

    i += 1
    print(f"{S} done")


print(Silh_scores)
Silh_scores_df = pd.DataFrame(
    np.transpose(Silh_scores),
    columns=["1.2", "1.4", "1.6", "1.8", "2.0", "2.2", "2.4", "2.6", "2.8", "3.0"],
)
Silh_scores_df.to_csv(f"silhouette_scores_{distance}.csv", index=False)
