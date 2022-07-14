import numpy as np


def generalized_distance_function(
    x_numeric, x_binary, x_categoric, distance_type, feature_counts
):

    if distance_type == "gower":
        out = gower_distance(x_numeric, x_categoric, x_binary, feature_counts)
    elif distance_type == "wishart":
        out = wishart_distance(x_numeric, x_categoric, x_binary, feature_counts)
    elif distance_type == "podani":
        out = podani_distance(x_numeric, x_categoric, x_binary, feature_counts)
    elif distance_type == "huang":
        out = huang_distance(x_numeric, x_categoric, x_binary)
    else:
        raise ValueError(
            "Not one of the supported distance method types: gower | wishart | podani | huang"
        )

    return out


def weighted_distance_matrix(
    per_feature_distances, weights, distance_type, feature_counts
):

    p = per_feature_distances.shape[0]
    n = per_feature_distances.shape[1]
    X = np.zeros((p, n, n))
    X[:] = per_feature_distances

    if distance_type == "gower" or distance_type == "huang":
        omega = 1.0
    elif distance_type == "podani" or distance_type == "wishart":
        omega = 0.5
    else:
        raise ValueError(
            "Not one of the supported distance method types: gower | wishart | podani | huang"
        )

    for k in range(p):  # Apply the feature weights to the distance matrix
        if k < feature_counts["Numeric"]:
            if distance_type == "gower":
                X[k, :, :] = X[k, :, :] * weights[k]
            else:
                X[k, :, :] = X[k, :, :] * weights[k] ** 2
                # Because Huang, Podani, and Wishart use a squared euclidean distance,
                # the weights must also be squared following the associative property of multiplication
        elif k < feature_counts["Numeric"] + feature_counts["Binary"]:
            X[k, :, :] = X[k, :, :] * weights[k]
        else:
            X[k, :, :] = X[k, :, :] * weights[k]

    out = sum(X) ** omega  # Sum along the k features and raise to power omega

    return out


### Heterogeneous data type distance functions ###


def gower_distance(x_numeric, x_categoric, x_binary, feature_counts):

    alpha = 1 / (sum(feature_counts.values()))
    beta = feature_counts["Binary"] / sum(feature_counts.values())
    gamma = feature_counts["Categoric"] / sum(feature_counts.values())

    out_numeric = alpha * range_weighted_manhattan(x_numeric)
    out_binary = beta * simple_matching(x_binary)
    out_categoric = gamma * simple_matching(x_categoric)
    out = np.concatenate(
        [x for x in [out_numeric, out_binary, out_categoric] if x.size > 0]
    )

    return out


def wishart_distance(x_numeric, x_categoric, x_binary, feature_counts):

    alpha = 1 / sum(feature_counts.values())
    beta = feature_counts["Binary"] / sum(feature_counts.values())
    gamma = feature_counts["Categoric"] / sum(feature_counts.values())

    out_numeric = alpha * variance_weighted_squared_euclidean(x_numeric)
    out_binary = beta * simple_matching(x_binary)
    out_categoric = gamma * simple_matching(x_categoric)
    out = np.concatenate(
        [x for x in [out_numeric, out_binary, out_categoric] if x.size > 0]
    )

    return out


def podani_distance(x_numeric, x_categoric, x_binary, feature_counts):

    alpha = 1.0
    beta = feature_counts["Binary"]
    gamma = feature_counts["Categoric"]

    out_numeric = alpha * squared_range_weighted_squared_euclidean(x_numeric)
    out_binary = beta * simple_matching(x_binary)
    out_categoric = gamma * simple_matching(x_categoric)
    out = np.concatenate(
        [x for x in [out_numeric, out_binary, out_categoric] if x.size > 0]
    )

    return out


def huang_distance(x_numeric, x_categoric, x_binary):

    alpha = 1.0
    beta = np.mean(np.std(x_numeric, axis=0, ddof=1))
    gamma = np.mean(np.std(x_numeric, axis=0, ddof=1))

    out_numeric = alpha * squared_euclidean(x_numeric)
    out_binary = beta * hamming(x_binary)
    out_categoric = gamma * hamming(x_categoric)
    out = np.concatenate(
        [x for x in [out_numeric, out_binary, out_categoric] if x.size > 0]
    )

    return out


### Homogeneous data type distance functions ###


def range_weighted_manhattan(X):

    if X.shape[1] == 0:
        return np.array([])

    out = np.zeros((X.shape[1], X.shape[0], X.shape[0]))
    numeric_ranges = np.zeros((X.shape[1]))
    for k in range(X.shape[1]):
        feature = X[:, k].astype(np.float32)
        numeric_ranges[k] = abs(np.nanmax(feature) - np.nanmin(feature))
    for i in range(len(X)):
        distance = np.abs(X[i, :] - X) / numeric_ranges
        out[:, i, :] = np.transpose(distance)
    return out


def squared_euclidean(X):

    if X.shape[1] == 0:
        return np.array([])

    out = np.zeros((X.shape[1], X.shape[0], X.shape[0]))

    for i in range(len(X)):
        distance = (X[i, :] - X) ** 2
        out[:, i, :] = np.transpose(distance)

    return out


def squared_range_weighted_squared_euclidean(X):

    if X.shape[1] == 0:
        return np.array([])

    out = np.zeros((X.shape[1], X.shape[0], X.shape[0]))
    numeric_ranges = np.zeros((X.shape[1]))
    for k in range(X.shape[1]):
        feature = X[:, k].astype(np.float32)
        numeric_ranges[k] = abs(np.nanmax(feature) - np.nanmin(feature))
    for i in range(len(X)):
        distance = ((X[i, :] - X) ** 2) / (numeric_ranges ** 2)
        out[:, i, :] = np.transpose(distance)
    return out


def variance_weighted_squared_euclidean(X):

    if X.shape[1] == 0:
        return np.array([])

    out = np.zeros((X.shape[1], X.shape[0], X.shape[0]))
    variance_weights = np.var(X, axis=0, ddof=1)
    for i in range(len(X)):
        distance = ((X[i, :] - X) ** 2) / variance_weights
        out[:, i, :] = np.transpose(distance)
    return out


def simple_matching(X):

    if X.shape[1] == 0:
        return np.array([])

    out = np.zeros((X.shape[1], X.shape[0], X.shape[0]))

    for i in range(len(X)):
        distance = (1.0 - 1.0 * (X[i, :] == X)) / X.shape[1]
        out[:, i, :] = np.transpose(distance)

    return out


def hamming(X):

    if X.shape[1] == 0:
        return np.array([])

    out = np.zeros((X.shape[1], X.shape[0], X.shape[0]))

    for i in range(len(X)):
        distance = 1.0 - 1.0 * (X[i, :] == X)
        out[:, i, :] = np.transpose(distance)

    return out
