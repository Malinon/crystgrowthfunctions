import unittest
import crystgrowthpoly.public as cryst
from sage.all import *


def general_testing_function(test_case_instance, symmetric_frame, file_path, starting_points, functions_triples):
    tes = cryst.read_tessellation_from_file(file_path)
    for i in range(len(starting_points)):
        growth_polynomials = tes.get_growth_polynomials_parallelogram(x0=starting_points[i], symmetric_frame=symmetric_frame)
        for j in range(3):
            test_case_instance.assertTrue(growth_polynomials[j].polynomials[0] == functions_triples[i][j])
class TestSimpleFrameGrowth(unittest.TestCase):

    def test_4444_symmetric(self):
        staring_points = ((0,0), (1/2, 1/3))
        var('n')
        functions_triples = ((n**2 + 2 * n + 1,  2*n**2 + 2*n,  n**2),
                            (n**2, 2*n**2 - 2*n, n**2 - 2*n + 1))
        general_testing_function(self, True, "./input/4444.txt", staring_points, functions_triples)
        # Asymetric case
        var("k l")
        functions_triples = ((k*l + k + l + 1, 2*k*l + k + l, k*l),
                           (k*l, 2*k*l - k - l, k*l - k - l + 1))
        general_testing_function(self, False, "./input/4444.txt", staring_points, functions_triples)

    def test_333333_symmetric(self):
        staring_points = ((0,0), (2/3, 1/2))
        var('n')
        functions_triples = ((2*n**2 + 2 * n + 1,  3*n**2,  n**2 - 2 * n + 1),
                            (2*n**2, 3*n**2 - 2*n, n**2 - 2*n + 1))
        general_testing_function(self, True, "./input/333333.txt", staring_points, functions_triples)
        # Asymetric case
        var("k l")
        functions_triples = ((2*k*l+k + l + 1, 3*k*l , k*l - k - l + 1),
                           (2*k*l, 3*k*l - k - l, k*l - k - l + 1))
        general_testing_function(self, False, "./input/333333.txt", staring_points, functions_triples)

    def test_pg_symmetric(self):
        staring_points = ((0,0), (1/4, 0))
        var('n')
        functions_triples = ((4*n**2 +4*n + 1,  6*n**2 + 3*n,  2*n**2 - n),
                            (4*n**2 + 2*n, 6*n**2 - 1, 2*n**2 - 2*n))
        general_testing_function(self, True, "./input/pg.txt", staring_points, functions_triples)
        # Asymetric case
        var("k l")
        functions_triples = ((4*k*l + 2*k + 2*l + 1, 6*k*l + k + 2*l, 2*k*l - k ),
                           (4*k*l + 2*k, 6*k*l + k - l - 1, 2*k*l - k - l))
        general_testing_function(self, False, "./input/pg.txt", staring_points, functions_triples)


if __name__ == '__main__':
    unittest.main()