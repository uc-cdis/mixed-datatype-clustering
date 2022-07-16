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
distance = "podani"
K = 5

hyperparams = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]

cluster_labels = np.zeros((X.shape[0], len(hyperparams)))
feature_weights = np.zeros((X.shape[1], len(hyperparams)))

i = 0
# for distance in distance_types:
for S in hyperparams:

    """
   spectral_cluster_labels, spectral_feature_weights, spectral_feature_order = clustering.spectral_kmedoids(
       X,
       distance_type=distance,
       k=K,
       method="pam",
       init="build",
       max_iter=100,
       random_state=None,
   )
   cluster_labels[:, K-2] = spectral_cluster_labels
   feature_weights[:, K-2] = spectral_feature_weights
   
   cluster_labels_df = pd.DataFrame(cluster_labels, columns=range(2, max_K + 1))
   cluster_labels_df.to_csv(f'spectral_cluster_labels_{distance}.csv', index=False)
   
   feature_weights_df = pd.DataFrame(feature_weights, columns=range(2, max_K + 1))
   feature_weights_df.to_csv(f'spectral_feature_weights_{distance}.csv', index=False)
   
   print(f'{distance} for k={K} done')
   
   """

    (
        sparse_cluster_labels,
        sparse_feature_weights,
        sparse_weight_difference,
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

    print(f"{S} done")

cluster_labels_df = pd.DataFrame(
    cluster_labels,
    columns=["1.2", "1.4", "1.6", "1.8", "2.0", "2.2", "2.4", "2.6", "2.8", "3.0"],
)
cluster_labels_df.to_csv(f"sparse_cluster_labels_{distance}.csv", index=False)

feature_weights_df = pd.DataFrame(
    feature_weights,
    columns=["1.2", "1.4", "1.6", "1.8", "2.0", "2.2", "2.4", "2.6", "2.8", "3.0"],
)
feature_weights_df.to_csv(f"sparse_feature_weights_{distance}.csv", index=False)