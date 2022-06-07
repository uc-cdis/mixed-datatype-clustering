# A Collection of k-medoid clustering and feature selection methods for mixed datatypes 

The first method is Sparse KMedoids Clustering which can use either of scikit-learn KMedoids 
clustering algorithms ('pam' or 'alternate'). This k-medoids method is based off of the sparse 
k-means clustering algorithm defined in this paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2930825/

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

#### clustering.sparse_kmedoids

    1) Calculates a generalized NxNxP featurewise distance matrix
    2) Calculate the weighted distance matrix using evenly weighted features
    3) Iteratively alternates between:
        i) Clustering 
        ii) Recalculating feature weights and weighted distance matrix 

    Returns the feature weights, cluster labels, the summed differences in feature weights on the final iteration, 
    and the order of the features in the feature weight vector. 

    The order of the features is recorded because regardless of the order the P features are given in the 
    original data matrix, the new data matrix is rearraged so that all of the numerical features are first, followed by
    the binary features and then the categorical features. The order of features is maintained amongst each relative
    feature type (i.e. numerical features maintain their order relative to each other).

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

    
