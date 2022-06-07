# Develop Sparse PAM Clustering by substituting PAM (k-medoids-based) clustering for k-means clustering in sparse k-means clustering algorithm  mentioned in this paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2930825/

## Some important notes on the functions and Sparse K-Medoids method

### Main Function: 

    Arguements:
    X -> The original data matrix, conintuous data must numeric
    distance_type -> The type of distance formula used. Currently supports gower, wishart, and podani distances
    k -> number of desired clusters
    s -> L1 norm upper bound of the feature weights (s >= ||w|| > 1.0) if s = None then performs tuning step first
    max_attempts -> The maximum number of attempts between clustering and feature selection
    method -> The method used in sklearn_extra.cluster.KMedoids for clustering ('pam' or 'alternate')
    init -> How the medoids are initialized in in sklearn_extra.cluster.KMedoids
    max_iter -> Max iterations used in sklearn_extra.cluster.KMedoids
    random_state -> If init = 'random' the random state seed used in sklearn_extra.cluster.KMedoids 

#### sparse_kmedoids

    1) Calculates a generalized NxNxP featurewise distance matrix
    2) Calculate the weighted distance matrix using evenly weighted features
    3) Iteratively alternates between:
        i) Clustering 
        ii) Recalculating feature weights and weighted distance matrix 

### Important Sub Functions:

#### 1) distfuncs.generalized_distance_function

    A master function for organizing the calculation of the NxNxP featureewise distance matrix

#### 2) distfuncs.weighted_distance_matrix

      Takes the NxNxP featurewise distance matrix and feature weights and returns the weighted distance matrix

#### 3) distfuncs.gower_distance

    Calculates the NxNxP gower distance matrix

#### 4) distfuncs.podani_distance

    Calculates the NxNxP podani distance matrix

#### 5) distfuncs.wishart_distance

    Calculates the NxNxP wishart distance matrix

#### 6) subfuncs.kmedoid_clusters

    Use the weighted distance matrix to perform kmedoid clustering and return new cluster labels

#### 7) subfuncs.between_medoid_sum_distances

    Calculates the difference between the total sums of distances and the sum of within cluster sums of distances

#### 8) subfuncs.update_weights

    Use the featurewise distance matrix and cluster labels to recalculate feature weights

    
