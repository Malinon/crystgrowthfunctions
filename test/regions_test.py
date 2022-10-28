import unittest
import crystgrowthpoly as cryst
import crystgrowthpoly.visualisation as cryst_visualisation
from sage.all import *


def general_testing_function(test_case_instance, file_path, parallelograms_number, edges_number,
                             vertices_number, color_numbers):
    tes = cryst.read_tessellation_from_file(file_path)
    regions = cryst_visualisation.find_regions(tes, symmetric_growth=True, full_plot=True)

    test_case_instance.assertTrue(len(regions.description_dict.keys()) == color_numbers,
                                  "Incorrect number of differnt growth functions")
    test_case_instance.assertEqual(parallelograms_number, len(regions.regions_lists[0]))
    test_case_instance.assertEqual(edges_number, len(regions.regions_lists[1]))
    test_case_instance.assertEqual(vertices_number, len(regions.regions_lists[2]))

class TestSimpleFrameGrowth(unittest.TestCase):

    def test_4444(self):
        general_testing_function(self, "./input/4444.txt", 1, 4, 4, 3)


    def test_333333(self):
        general_testing_function(self, "./input/333333.txt", 4, 12, 9, 7)


if __name__ == '__main__':
    unittest.main()