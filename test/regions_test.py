import unittest
import crystgrowthpoly as cryst
import crystgrowthpoly.visualisation as cryst_visualisation
from sage.all import *


def general_testing_function(test_case_instance, file_path, parallelograms_number, edges_number,
                             vertices_number, color_numbers):
    tes = cryst.read_tessellation_from_file(file_path)
    regions = cryst_visualisation.find_regions(tes, symmetric_frame=True, full_plot=True)

    test_case_instance.assertTrue(regions.parallelograms_plot != None)
    test_case_instance.assertTrue(regions.edges_plot != None)
    test_case_instance.assertTrue(regions.corners_plot != None)
    test_case_instance.assertTrue(regions.regions_description != None)
    par_counter = 0
    edge_counter = 0
    vertices_counter = 0
    test_case_instance.assertTrue(len(regions.regions_description.keys()) == color_numbers,
                                  "Incorrect number of differnt growth functions")
    for key in regions.regions_description.keys():
        for region in regions.regions_description[key]:
            if isinstance(region, cryst_visualisation.OpenRectangle):
                par_counter = par_counter + 1
            elif isinstance(region, cryst_visualisation.OpenLine):
                edge_counter = edge_counter + 1
            else:
                test_case_instance.assertTrue(isinstance(region, cryst_visualisation.SpecialPoint), "Incorrect type")
                vertices_counter = vertices_counter + 1
    test_case_instance.assertEqual(parallelograms_number, par_counter)
    test_case_instance.assertEqual(edges_number, edge_counter)
    test_case_instance.assertEqual(vertices_number, vertices_counter)

class TestSimpleFrameGrowth(unittest.TestCase):

    def test_4444(self):
        general_testing_function(self, "./input/4444.txt", 1, 4, 4, 3)


    def test_333333(self):
        general_testing_function(self, "./input/333333.txt", 4, 12, 9, 7)


if __name__ == '__main__':
    unittest.main()