import unittest
import crystgrowthpoly.public as cryst
from sage.all import *



def general_testing_function(test_case_instance, symmetric_frame, file_path, functions_triple):
    tes = cryst.read_tessellation_from_file(file_path)
    growth_polynomials = tes.get_growth_polynomials(symmetric_frame)
    for j in range(3):
        test_case_instance.assertTrue(growth_polynomials[j].polynomials[0] == functions_triple[j])

class TestTopologicalGrowth(unittest.TestCase):

    def test_4444(self):
        var('n')
        functions_triple = (n**2 + 2 * n + 1,  2*n**2 + 2*n,  n**2)
        general_testing_function(self, True, "./input/4444.txt", functions_triple)
        # Asymetric case
        var("n_1 n_2")
        functions_triple = (n_1*n_2 + n_1 + n_2 + 1, 2*n_1*n_2 + n_1 + n_2, n_1*n_2)
        general_testing_function(self, False, "./input/4444.txt", functions_triple)

    def test_333333(self):
        var('n')
        functions_triple = (2 * n**2 + 4 * n,  3*n**2 + 4*n -1,  n**2)
        general_testing_function(self, True, "./input/333333.txt", functions_triple)
        # Asymetric case
        var("n_1 n_2")
        functions_triple = (2*n_1*n_2 + 2*n_1 + 2*n_2, 3*n_1*n_2  + 2*n_1  + 2*n_2 - 1, n_1*n_2)
        general_testing_function(self, False, "./input/333333.txt", functions_triple)

    def test_p2(self):
        var('n')
        functions_triple = (2*n**2 + 3*n + 1, 4*n**2 + 3*n, 2*n**2)
        general_testing_function(self, True, "./input/p2.txt", functions_triple)
        # Asymetric case
        var("n_1 n_2")
        functions_triple = (2*n_1*n_2 + 2*n_1 + n_2 + 1, 4*n_1*n_2 + 2*n_1 + n_2, 2*n_1*n_2)
        general_testing_function(self, False, "./input/p2.txt", functions_triple)

    def test_pg(self):
        var('n')
        functions_triple = (4*n**2 + 6*n, 6*n**2 + 6*n - 1, 2 * n**2)
        general_testing_function(self, True, "./input/pg.txt", functions_triple)
        # Asymetric case
        var("n_1 n_2")
        functions_triple = (4*n_1*n_2 + 4*n_1 + 2*n_2, 6*n_1*n_2 + 4*n_1 + 2*n_2 - 1, 2 * n_1 * n_2)
        general_testing_function(self, False, "./input/pg.txt", functions_triple)


if __name__ == '__main__':
    unittest.main()