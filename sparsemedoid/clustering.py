from sparsemedoid.distfuncs import (
    weighted_distance_matrix,
    generalized_distance_function,
)
from sparsemedoid.subfuncs import (
    sort_datatypes,
    kmedoid_clusters,
    update_weights,
    spectral_feature_selection,
)

import numpy as np


def sparse_kmedoids(
    X, distance_type, k, s, max_attempts, method, init, max_iter=None, random_state=None
):

    if init == "random" and random_state is None:
        raise ValueError("Using `random` with no `random_state` seed")
    if s <= 1.0:
        raise ValueError("Hyperparameter s must be > 1.0")

    n, p = X.shape

    x_numeric, x_binary, x_categoric, feature_order = sort_datatypes(X, p)
    feature_counts = {
        "Numeric": x_numeric.shape[1],
        "Binary": x_binary.shape[1],
        "Categoric": x_categoric.shape[1],
    }

    per_feature_distances = generalized_distance_function(
        x_numeric, x_binary, x_categoric, distance_type, feature_counts
    )

    weights = np.ones(p) * (1 / np.sqrt(p))
    weighted_distances = weighted_distance_matrix(
        per_feature_distances, weights, distance_type, feature_counts
    )
    cluster_labels = kmedoid_clusters(
        weighted_distances, k, method, init, max_iter, random_state
    )

    attempts = 0
    weight_difference = 1.0

    while weight_difference > 1e-4 and attempts <= max_attempts:

        previous_weights = weights

        if attempts == 0:
            weights = update_weights(per_feature_distances, cluster_labels, s)
        else:
            weighted_distances = weighted_distance_matrix(
                per_feature_distances, weights, distance_type, feature_counts
            )
            cluster_labels = kmedoid_clusters(
                weighted_distances, k, method, init, max_iter, random_state
            )
            weights = update_weights(per_feature_distances, cluster_labels, s)

        attempts += 1
        weight_difference = np.sum(np.abs(weights - previous_weights)) / np.sum(
            np.abs(previous_weights)
        )

    return cluster_labels, weights, feature_order, weighted_distances


def spectral_kmedoids(
    X, distance_type, k, method, init, max_iter=None, random_state=None
):
    if init == "random" and random_state is None:
        raise ValueError("Using `random` with no `random_state` seed")

    n, p = X.shape

    x_numeric, x_binary, x_categoric, feature_order = sort_datatypes(X, p)
    feature_counts = {
        "Numeric": x_numeric.shape[1],
        "Binary": x_binary.shape[1],
        "Categoric": x_categoric.shape[1],
    }

    per_feature_distances = generalized_distance_function(
        x_numeric, x_binary, x_categoric, distance_type, feature_counts
    )

    unordered_weights = spectral_feature_selection(
        per_feature_distances, k, distance_type, feature_counts
    )

    weights = np.zeros(p)
    for i in range(p):
        weights[i] = unordered_weights[str(i)]

    weighted_distances = weighted_distance_matrix(
        per_feature_distances, weights, distance_type, feature_counts
    )

    cluster_labels = kmedoid_clusters(
        weighted_distances, k, method, init, max_iter, random_state
    )

    return cluster_labels, weights, feature_order
