import crystgrowthpoly.private as cryst_private
import crystgrowthpoly.visualisation as cryst_visualisation
import crystgrowthpoly.growth_function as gf
import re
from sage.all import *

"""
Function finding growth polynomials

@param cells Sorted (by cells dimmension) list of cells list
@param translation_vectors  Vectors generating tessellation
"""
def get_topological_growth_polynomials(cells, translation_vectors, symmetric_frame=True, normalized=True):
    dim = len(cells) - 1
    if symmetric_frame:
        args = tuple(tuple(i for j in range(dim)) for i in range(1, dim + 2))
        variables_num = 1
    else:
        # TODO: Implement generating arguments for higher dimmensions.
        if dim == 2:
            args = ((1,1), (1,2), (2,1), (2,2))
        elif dim == 3:
            args = ((1,1,1), (1,1,2), (1,2,1), (2,1,1), (2,2,1), (2,1,2), (1,2,2), (2,2,2))
        else:
            raise NotImplementedError("It is not implemented yet.")
        variables_num = dim
    polynomial_k_cells_coeff_lists = list()
    if not normalized:
        args = tuple(tuple(2 * val for val in arg) for arg in args)
    for i in range(len(cells) - 2):
        gen_val = lambda arg: cryst_private.get_k_cells_num(arg, cells[i], translation_vectors)
        polynomial_k_cells_coeff_lists.append(cryst_private.find_poly(gen_val, args, symmetric_frame)[0])
    # Calculate coefficients of growth function related to cells of highest dimmensions.
    highest_dim_function = list(0 for i in range(len(args)))
    highest_dim_function[0] = len(cells[-1])
    # Calculate coefficients of growth functions related to dim - 1 cells based on Euler Characteristic
    polynomial_k_cells_coeff_lists.append(
        cryst_private.calculate_coefficients_based_ne_euler_charcteristics(
            polynomial_k_cells_coeff_lists, highest_dim_function, dim))
    polynomial_k_cells_coeff_lists.append(highest_dim_function )
    return tuple(gf.growth_function((polynomial_k_cells_coeff_lists[i],), variables_num, i, normalized=normalized) for i in range(dim + 1))


"""
Function finding growth polynomials

@param points 0-cells of repeating figure
@param edges 1-cells of repeating figure
@param faces 2-cells of r
@param v1  tesselation vector
@param v2  tesselation vector
@param x0  
@param frame_scale_v1  frame_scale_v1 * v1 is used to define parallelogram frame
@param frame_scale_v2  frame_scale_v2 * v2 is used to define parallelogram frame
"""
def get_crystalographic_growth_functions(points, edges, faces, v1, v2, x0, frame_scale_v1 = 1, frame_scale_v2 = 1, in_one_variable=True):
    try:
        rational_scale_v1 = Rational(frame_scale_v1)
        rational_scale_v2 = Rational(frame_scale_v2)
    except TypeError:
        raise TypeError("Frame scale need to be rational number")
    # Find point equivalent to x0 near provided figure
    params = column_matrix([v1, v2]).solve_right(vector(cryst_private.subtract_vectors(points[0], x0)))
    x0_equivalent = cryst_private.translate_vector(cryst_private.translate_vector(x0, cryst_private.multiply_vector(v1, floor(params[0]))),
                                                   cryst_private.multiply_vector(v2, floor(params[1])))
    v1_frame = cryst_private.multiply_vector(v1, frame_scale_v1)
    v2_frame = cryst_private.multiply_vector(v2, frame_scale_v2)
    limits_extenders  = cryst_private.get_limits_extenders(x0_equivalent, points)
    polynomials_0_cells = []
    polynomials_1_cells = []
    polynomials_2_cells = []
    if in_one_variable:
        N = lcm(rational_scale_v1.denominator(), rational_scale_v2.denominator())
        # Prepare functions calculating numbers of cells
        gen_num_of_0_cells = lambda arg: cryst_private.get_k_cells_num_parallelogram(
            arg[0], arg[0], points, v1, v2, 0, x0_equivalent, ceil(frame_scale_v1), ceil(frame_scale_v2), v1_frame, v2_frame, limits_extenders)
        gen_num_of_1_cells = lambda arg: cryst_private.get_k_cells_num_parallelogram(
            arg[0], arg[0], edges, v1, v2, 1, x0_equivalent, ceil(frame_scale_v1), ceil(frame_scale_v2), v1_frame, v2_frame, limits_extenders)
        gen_num_of_2_cells = lambda arg: cryst_private.get_k_cells_num_parallelogram(
            arg[0], arg[0], faces, v1, v2, 2, x0_equivalent,ceil(frame_scale_v1), ceil(frame_scale_v2), v1_frame, v2_frame, limits_extenders)
        for s in range(N):
            args = tuple((N*i + s,) for i in range(2,5) )
            # Find growth polynomials
            polynomials_0_cells.append(cryst_private.find_poly(gen_num_of_0_cells, args, symmetric_frame=in_one_variable)[0])
            polynomials_1_cells.append(cryst_private.find_poly(gen_num_of_1_cells, args, symmetric_frame=in_one_variable)[0])
            polynomials_2_cells.append(cryst_private.find_poly(gen_num_of_2_cells, args, symmetric_frame=in_one_variable)[0])
        return (gf.growth_function(polynomials_0_cells, 1, 0, N), gf.growth_function(polynomials_1_cells, 1, 1, N), gf.growth_function(polynomials_2_cells, 1, 2, N))
    else:
        K = rational_scale_v1.denominator()
        L = rational_scale_v2.denominator()
        gen_num_of_0_cells = lambda arg: cryst_private.get_k_cells_num_parallelogram(
            arg[0], arg[1], points, v1, v2, 0, x0_equivalent, ceil(frame_scale_v1), ceil(frame_scale_v2), v1_frame, v2_frame, limits_extenders)
        gen_num_of_1_cells = lambda arg: cryst_private.get_k_cells_num_parallelogram(
            arg[0], arg[1], edges, v1, v2, 1, x0_equivalent, ceil(frame_scale_v1), ceil(frame_scale_v2), v1_frame, v2_frame, limits_extenders)
        gen_num_of_2_cells = lambda arg: cryst_private.get_k_cells_num_parallelogram(
            arg[0], arg[1], faces, v1, v2, 2, x0_equivalent,ceil(frame_scale_v1), ceil(frame_scale_v2), v1_frame, v2_frame, limits_extenders)
        for s_k in range(K):
            for s_l in range(L):
                args = tuple( (s_k + arg[0] * K, s_l + arg[1] * L)  for arg in ((2,2), (2,3), (3,2), (3,3)))
                polynomials_0_cells.append(cryst_private.find_poly(gen_num_of_0_cells, args, in_one_variable)[0])
                polynomials_1_cells.append(cryst_private.find_poly(gen_num_of_1_cells, args, in_one_variable)[0])
                polynomials_2_cells.append(cryst_private.find_poly(gen_num_of_2_cells, args, in_one_variable)[0])
        return (gf.growth_function(polynomials_0_cells, 2, 0, (K, L)), gf.growth_function(polynomials_1_cells, 2, 1, (K, L)), gf.growth_function(polynomials_2_cells, 2, 2, (K, L)))

class Polygon:
    def __init__(self, cells):
        self.cells = list()
        self.cells.append(cells[0])
        for i in range(1, len(cells)):
            self.cells.append(tuple(cryst_private.sort_points_in_cell(cell) for cell in cells[i]))

class Tesselation:
    __FILE_PREFIXES = ("domains_parallelogram", "domains_lines_and_points")
    def __init__(self, polygon, cartesian_tesselation_vectors=None):
        self.polygon = polygon
        self.translation_vectors =tuple(tuple(row) for row in matrix.identity(len(self.polygon.cells) - 1))
        self.dim = len(self.polygon.cells) - 1
        # TODO: Remove after implementing finding crystallographic growth functions for higher dimmensions
        self.v1 = (1,0)
        self.v2 = (0,1)
        self.cartesian_vectors = cartesian_tesselation_vectors
    def get_growth_polynomials_parallelogram(self, scale_v1 = 1, scale_v2 = 1, x0 = (0,0), symmetric_frame=True):
        return get_crystalographic_growth_functions(self.polygon.cells[0], self.polygon.cells[1], self.polygon.cells[2], self.v1, self.v2, x0,
                                                    frame_scale_v1 = scale_v1, frame_scale_v2 = scale_v2,
                                                   in_one_variable=symmetric_frame)
    def get_growth_polynomials(self, symmetric_frame=True, normalized = True):
        return get_topological_growth_polynomials(self.polygon.cells, self.translation_vectors, symmetric_frame, normalized)
    def plot_edges(self, repetition_of_polygon):
        if repetition_of_polygon == None:
            repetition_of_polygon = ((0,0), (0,0))
        coord_matrix = matrix(self.cartesian_vectors)
        get_edge_in_cartesian_coordinates = lambda edge: (vector(edge[0]) * coord_matrix, vector(edge[1]) * coord_matrix)
        G = Graphics()
        for i in range(repetition_of_polygon[0][0],  repetition_of_polygon[0][1] + 1):
            for j in range(repetition_of_polygon[1][0],  repetition_of_polygon[1][1] + 1):
                G = G + sum( line(get_edge_in_cartesian_coordinates(cryst_private.translate_face(e, (i,j))), color="black" )
                            for e in self.polygon.cells[1])
        return G
    def plot_domains(self, symmetric_growth = True, full_plot=False, description = True, export_format=None, file_name_sufix = "",
                     xmin=None, xmax=None, ymin=None, ymax=None, repetition_of_unit_cells = None, repetition_of_polygon = None, fit_image=False,
                     return_object=False, colors_lists=(None, None)):
        regions_desc = cryst_visualisation.find_regions(self, symmetric_growth, full_plot=full_plot) # repetition_of_unit_cells=repetition_of_unit_cells)
        if return_object:
            return regions_desc
        plots = regions_desc.get_plots(repetition_of_unit_cells=repetition_of_unit_cells, repetition_of_polygon=repetition_of_polygon,
                                       colors_lists=colors_lists)
        for i in range(len(plots)):
            plots[i].show(aspect_ratio=1, axes=False, xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax)
            if export_format != None:
                plots[i].save(__FILE_PREFIXES[i] + file_name_sufix + "." + export_format,
                              xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, axes=False)

        if description:
            regions_desc.describe()


def read_tessellation_from_file(file_path, cartesian_vectors_included=True, dim=2, crystallographic_coordinates=True):
    '''Reads tessellation from file

Input File Format:
    The first line shows the number of cells in each dimension sorted in ascending order with respect to the cell dimension.
    The following lines encode the coordinates of the tessellation vertices.
    Coordinates can be given in crystallographic or Cartesian coordinates.
    When describing vertices, the SageMath syntax can be used.
    Higher dimensional cells are described after the vertices.
    Each line encodes one cell.
    Higher dimensional cells are encoded with vertex indices.
    The last two lines contain information about translation vectors

Example input file:
    6 6 1
        0     0
        1/3  -1/3
        1     0
        4/3   2/3
        1     1
        1/3   2/3
    1 2
    2 3
    3 4
    4 5
    5 6
    6 1
    1 2 3 4 5 6
    sqrt(3)*2 0
    -sqrt(3) 3
    '''
    input_file = open(file_path, 'r')
    lines = input_file.readlines()
    input_file.close()
    cells_nums = tuple(int(s) for s in re.findall(r"\d+", lines[0]))
    cells = list()
    cartesian_vectors = None
    if cartesian_vectors_included:
        cartesian_vectors = tuple(cryst_private.string_to_point(lines[i-dim]) for i in range(dim))
    if crystallographic_coordinates:
        points = tuple(cryst_private.string_to_point(lines[i])
                  for i in range(1, cells_nums[0] + 1))
    else:
        converter = lambda p: tuple(vector(p) * ~matrix(cartesian_vectors))
        points = tuple(converter(cryst_private.string_to_point(lines[i]))
                  for i in range(1, cells_nums[0] + 1))
    cells.append(points)
    acc_cells_number = cells_nums[0] + 1
    for i in range(1, dim + 1):
        cells.append(tuple(cryst_private.string_to_cell(lines[j], points)
                  for j in range(acc_cells_number, acc_cells_number + cells_nums[i])))
        acc_cells_number = acc_cells_number + cells_nums[i]
    return Tesselation(Polygon(cells),cartesian_vectors)
