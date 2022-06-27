import unittest
import numpy as np
from sparsemedoid import distfuncs, subfuncs


class test_subfuncs(unittest.TestCase):
    def test_between_medoid_sum_distances(self):

        X_diss = (
            np.array(
                [
                    [[0.0, 1.0, 0.25], [1.0, 0.0, 0.75], [0.25, 0.75, 0.0]],
                    [[0.0, 0.5, 1.0], [0.5, 0.0, 0.5], [1.0, 0.5, 0.0]],
                    [[0.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
                ]
            )
            / 3.0
        )

        cls1 = np.array([1, 0, 1])
        actual1 = subfuncs.between_medoid_sum_distances(X_diss, cls1)[0]
        expected1 = np.array([13.0 / 36.0, 1.0 / 9.0, 4.0 / 9.0])

        cls2 = np.array([1, 1, 0])
        actual2 = subfuncs.between_medoid_sum_distances(X_diss, cls2)[0]
        expected2 = np.array([1.0 / 9.0, 5.0 / 18.0, 1.0 / 9.0])

        np.testing.assert_almost_equal(expected1, actual1)
        np.testing.assert_almost_equal(expected2, actual2)

    def test_positive_part(self):

        x1 = np.array([1, 0.5, 0.002, 1.04, 0])
        x2 = np.array([0.3, 2.4, -0.4, 0.005, 0])
        x3 = np.array([-1, -0.3, -1.6, -0.0004, -223])

        expected1 = np.array([1, 0.5, 0.002, 1.04, 0])
        expected2 = np.array([0.3, 2.4, 0, 0.005, 0])
        expected3 = np.array([0, 0, 0, 0, 0])

        actual1 = subfuncs.positive_part(x1)
        actual2 = subfuncs.positive_part(x2)
        actual3 = subfuncs.positive_part(x3)

        np.testing.assert_almost_equal(expected1, actual1)
        np.testing.assert_almost_equal(expected2, actual2)
        np.testing.assert_almost_equal(expected3, actual3)

    def test_scale_weights(self):

        x1 = np.array([1, 4, 3, 5])
        x2 = np.array([0, 0.2, 1.5, 3])

        expected1 = np.array(
            [1.0 / np.sqrt(51), 4 / np.sqrt(51), 3 / np.sqrt(51), 5 / np.sqrt(51)]
        )
        expected2 = np.array(
            [0, 0.2 / np.sqrt(11.29), 1.5 / np.sqrt(11.29), 3 / np.sqrt(11.29)]
        )

        actual1 = subfuncs.scale_weights(x1)
        actual2 = subfuncs.scale_weights(x2)

        np.testing.assert_almost_equal(expected1, actual1)
        np.testing.assert_almost_equal(expected2, actual2)

    def test_soft_thresholding(self):

        x1 = np.array([1, 2, -4, 0, 0.3, 1.4])
        d1 = 0.8

        x2 = np.array([-1, 3.2, -0.23, 10, 1.4, -1.9])
        d2 = 1.3

        actual1 = subfuncs.soft_thresholding(x1, d1)
        expected1 = np.array([0.2, 1.2, -3.2, 0, 0, 0.6])

        actual2 = subfuncs.soft_thresholding(x2, d2)
        expected2 = np.array([0, 1.9, 0, 8.7, 0.1, -0.6])

        np.testing.assert_almost_equal(expected1, actual1)
        np.testing.assert_almost_equal(expected2, actual2)

    def test_binary_search(self):

        a1 = np.array([0, 0, 0, 0])
        a2 = np.array([1, 2, 3, 4])
        a3 = np.array([0.0, 0.5, 0.5, 2 / 3, 1.0])

        s1 = 5.0
        s2 = 5.0
        s3 = 1.5

        expected1 = 0.0
        expected2 = 0.0
        expected3 = 0.4352068

        actual1 = subfuncs.binary_search(a1, s1)
        actual2 = subfuncs.binary_search(a2, s2)
        actual3 = subfuncs.binary_search(a3, s3)

        self.assertEqual(expected1, actual1)
        self.assertEqual(expected2, actual2)

        self.assertAlmostEqual(expected3, actual3)

    def test_simple_matching(self):

        X1 = np.array([[]])
        X2 = np.array([["A", "1"], ["B", "1"], ["A", "2"]])

        expected1 = np.array([])
        actual1 = distfuncs.simple_matching(X1)

        expected2 = np.array(
            [
                [[0, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0]],
                [[0, 0, 0.5], [0, 0, 0.5], [0.5, 0.5, 0]],
            ]
        )
        actual2 = distfuncs.simple_matching(X2)

        np.testing.assert_almost_equal(expected1, actual1)
        np.testing.assert_almost_equal(expected2, actual2)

    def test_hamming(self):

        X1 = np.array([[]])
        X2 = np.array([["A", "1"], ["B", "1"], ["A", "2"]])

        expected1 = np.array([])
        actual1 = distfuncs.hamming(X1)

        expected2 = np.array(
            [[[0, 1, 0], [1, 0, 1], [0, 1, 0]], [[0, 0, 1], [0, 0, 1], [1, 1, 0]]]
        )
        actual2 = distfuncs.hamming(X2)

        np.testing.assert_almost_equal(expected1, actual1)
        np.testing.assert_almost_equal(expected2, actual2)

    def test_range_weighted_manhattan(self):
        X1 = np.array([[]])
        X2 = np.array([[0.2, 0.0], [1.0, 0.7], [0.0, 1.0]])

        expected1 = np.array([])
        actual1 = distfuncs.range_weighted_manhattan(X1)

        expected2 = np.array(
            [
                [[0, 0.8, 0.2], [0.8, 0, 1.0], [0.2, 1.0, 0]],
                [[0, 0.7, 1.0], [0.7, 0, 0.3], [1.0, 0.3, 0]],
            ]
        )
        actual2 = distfuncs.range_weighted_manhattan(X2)

        np.testing.assert_almost_equal(expected1, actual1)
        np.testing.assert_almost_equal(expected2, actual2)

    def test_squared_euclidean(self):
        X1 = np.array([[]])
        X2 = np.array([[0.2, 0.0], [1.0, 0.7], [0.0, 1.0]])

        expected1 = np.array([])
        actual1 = distfuncs.squared_euclidean(X1)

        expected2 = np.array(
            [
                [[0, 0.64, 0.04], [0.64, 0, 1.0], [0.04, 1.0, 0]],
                [[0, 0.49, 1.0], [0.49, 0, 0.09], [1.0, 0.09, 0]],
            ]
        )

        actual2 = distfuncs.squared_euclidean(X2)

        np.testing.assert_almost_equal(expected1, actual1)
        np.testing.assert_almost_equal(expected2, actual2)

    def test_squared_range_weighted_squared_euclidean(self):
        X1 = np.array([[]])
        X2 = np.array([[0.2, 0.0], [1.0, 0.7], [0.0, 1.0]])

        expected1 = np.array([])
        actual1 = distfuncs.squared_range_weighted_squared_euclidean(X1)

        expected2 = np.array(
            [
                [[0, 0.64, 0.04], [0.64, 0, 1.0], [0.04, 1.0, 0]],
                [[0, 0.49, 1.0], [0.49, 0, 0.09], [1.0, 0.09, 0]],
            ]
        )
        actual2 = distfuncs.squared_range_weighted_squared_euclidean(X2)

        np.testing.assert_almost_equal(expected1, actual1)
        np.testing.assert_almost_equal(expected2, actual2)

    def test_variance_weighted_squared_euclidean(self):
        X1 = np.array([[]])
        X2 = np.array([[0.2, 0.0], [1.0, 0.7], [0.0, 1.0]])

        expected1 = np.array([])
        actual1 = distfuncs.variance_weighted_squared_euclidean(X1)

        expected2 = np.array(
            [
                [
                    [0, 2.2857143, 0.1428571],
                    [2.2857143, 0, 3.5714286],
                    [0.1428571, 3.5714286, 0],
                ],
                [
                    [0, 1.8607595, 3.7974684],
                    [1.8607595, 0, 0.34177216],
                    [3.7974684, 0.34177216, 0],
                ],
            ]
        )
        actual2 = distfuncs.variance_weighted_squared_euclidean(X2)

        np.testing.assert_almost_equal(expected1, actual1)
        np.testing.assert_almost_equal(expected2, actual2)


if __name__ == "__main__":
    unittest.main()
