import numpy as np
from sklearn_extra.cluster import KMedoids


def sort_datatypes(X, p):

    # FEATURE LABELS #
    # 0 = Numerical Data
    # 1 = Categorical Data
    # 2 = Binary Data

    feature_labels = np.zeros(p)
    for i in range(p):
        if type(X[0, i]) is str:
            if len(np.unique(X[:, i])) > 2:
                feature_labels[i] = 1
            else:
                feature_labels[i] = 2

    # Sort and group all data into 3 new dataframes. One dataframe for each datatype (numeric, categoric, and binary)

    numeric_columns = list()
    categoric_columnms = list()
    binary_columns = list()
    for i in range(p):
        if feature_labels[i] == 0:
            numeric_columns.append(i)
        if feature_labels[i] == 1:
            categoric_columnms.append(i)
        if feature_labels[i] == 2:
            binary_columns.append(i)

    x_numeric = np.zeros((len(X), len(numeric_columns)))
    x_binary = np.zeros((len(X), len(binary_columns))).astype(str)
    x_categoric = np.zeros((len(X), len(categoric_columnms))).astype(str)

    x_numeric[:] = X[:, numeric_columns]
    x_binary[:] = X[:, binary_columns]
    x_categoric[:] = X[:, categoric_columnms]

    feature_order = {
        "Numerical Features": numeric_columns,
        "Binary Features": binary_columns,
        "Categorical Features": categoric_columnms,
    }

    return x_numeric, x_binary, x_categoric, feature_order


def kmedoid_clusters(weighted_distances, k, method, init, max_iter, random_state):

    kmds = KMedoids(
        n_clusters=k,
        init=init,
        max_iter=max_iter,
        metric="precomputed",
        method=method,
        random_state=random_state,
    ).fit(weighted_distances)
    cluster_labels = kmds.labels_

    return cluster_labels


def update_weights(per_feature_distances, cluster_labels, s):

    p = per_feature_distances.shape[0]
    n = per_feature_distances.shape[1]
    X = np.zeros((p, n, n))
    X[:] = per_feature_distances

    bmsd = between_medoid_sum_distances(X, cluster_labels)[0]
    bmsd = positive_part(bmsd)
    delta = binary_search(bmsd, s)
    unscaled_weights = soft_thresholding(bmsd, delta)
    scaled_weights = scale_weights(unscaled_weights)

    return scaled_weights


def between_medoid_sum_distances(per_feature_distances, cluster_labels):

    p = per_feature_distances.shape[0]
    n = per_feature_distances.shape[1]
    X = np.zeros((p, n, n))
    X[:] = per_feature_distances

    tot_dis = np.sum(np.sum(X, axis=1), axis=1) / n
    cls_dis = np.zeros(p)

    for i in np.unique(cluster_labels):
        mask = cluster_labels == i
        n_k = np.sum(mask)
        if n_k > 1:  # If n_k = 1 then the within medoid distance is 0
            cls_dis += np.sum(np.sum(X[:, mask, :], axis=1)[:, mask], axis=1) / n_k

    a = tot_dis - cls_dis
    bmsd = np.sum(a)

    return a, bmsd


def soft_thresholding(bmsd, delta):

    threshold = np.sign(bmsd) * np.maximum(0, np.abs(bmsd) - delta)
    return threshold


def binary_search(bmsd, s):

    l2n_argu = np.linalg.norm(bmsd)  # 6
    if l2n_argu == 0 or np.sum(np.abs(bmsd / l2n_argu)) <= s:
        return 0
    lam1 = 0
    lam2 = np.max(np.abs(bmsd)) - 1e-5  # 0.99999
    i = 1

    while i <= 15 and (lam2 - lam1) > 1e-4:
        su = soft_thresholding(bmsd, (lam1 + lam2) / 2.0)
        if np.sum(np.abs(su / np.linalg.norm(su))) < s:
            lam2 = (lam1 + lam2) / 2.0
        else:
            lam1 = (lam1 + lam2) / 2.0
        i += 1

    return (lam1 + lam2) / 2.0


def positive_part(bmsd):

    positive_bmsd = [x if x >= 0 else 0 for x in bmsd]  # Take positive part of a
    return positive_bmsd


def scale_weights(unscaled_weights):

    return unscaled_weights / np.linalg.norm(unscaled_weights)
