""" Functions"""
import crystgrowthpoly.private as cryst_private
import crystgrowthpoly.visualisation as cryst_visualisation
import crystgrowthpoly.growth_function as gf
import re
from sage.all import *

## @param cells Sorted (by cells dimmension) list of cells list
## @param translation_vectors  Vectors generating tessellation
## @param symmetric_growth True, if growth is the same in each direction. False otherwise.
## @param normalized Normalization can be described as dividing argument by 2
## @return Triple of growth functions @ref growth_function . Functions are sorted ascending by dimmensions of corresponding cells.
def get_topological_growth_polynomials(cells, translation_vectors, symmetric_growth=True, normalized=True):
"""
Function finding topological growth functions
"""
    dim = len(cells) - 1
    if symmetric_growth:
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
    # Calculate coefficients of growth function for cells with dimmension from 0 to dim - 2
    for i in range(len(cells) - 2):
        gen_val = lambda arg: cryst_private.get_k_cells_num(arg, cells[i], translation_vectors)
        polynomial_k_cells_coeff_lists.append(cryst_private.find_poly(gen_val, args, symmetric_growth)[0])
    # Calculate coefficients of growth function related to cells of highest dimmensions.
    highest_dim_function = list(0 for i in range(len(args)))
    highest_dim_function[0] = len(cells[-1])
    # Calculate coefficients of growth functions related to dim - 1 cells based on Euler Characteristic
    polynomial_k_cells_coeff_lists.append(
        cryst_private.calculate_coefficients_based_ne_euler_charcteristics(
            polynomial_k_cells_coeff_lists, highest_dim_function, dim))
    polynomial_k_cells_coeff_lists.append(highest_dim_function)
    return tuple(gf.growth_function((polynomial_k_cells_coeff_lists[i],), variables_num, i, normalized=normalized) for i in range(dim + 1))


## @param points 0-cells of repeating figure
## @param edges 1-cells of repeating figure
## @param faces 2-cells of r
## @param v1  Tessellation vector
## @param v2  Tessellation vector
## @param x0  
## @param frame_scale_v1  frame_scale_v1 * v1 is used to define parallelogram frame
## @param frame_scale_v2  frame_scale_v2 * v2 is used to define parallelogram frame
## @param symmetric_growth True, if growth is the same in each direction. False otherwise.
def get_crystallographic_growth_functions(points, edges, faces, v1, v2, x0, frame_scale_v1 = 1, frame_scale_v2 = 1, symmetric_growth=True):
"""
Function finding crystallographic_growth_functions
"""
    try:
        rational_scale_v1 = Rational(frame_scale_v1)
        rational_scale_v2 = Rational(frame_scale_v2)
    except TypeError:
        raise TypeError("Frame scale need to be rational number")
    # Find point equivalent to x0 near provided figure
    params = column_matrix([v1, v2]).solve_right(vector(cryst_private.subtract_vectors(points[0], x0)))
    x0_equivalent = cryst_private.translate_vector(cryst_private.translate_vector(x0, cryst_private.multiply_vector(v1, floor(params[0]))),
                                                   cryst_private.multiply_vector(v2, floor(params[1])))
    # Scaled vectors
    v1_frame = cryst_private.multiply_vector(v1, frame_scale_v1)
    v2_frame = cryst_private.multiply_vector(v2, frame_scale_v2)
    limits_extenders  = cryst_private.get_limits_extenders(x0_equivalent, points)
    polynomials_0_cells = []
    polynomials_1_cells = []
    polynomials_2_cells = []
    if symmetric_growth:
        N = lcm(rational_scale_v1.denominator(), rational_scale_v2.denominator())
        # Prepare functions calculating numbers of cells
        gen_num_of_0_cells = lambda arg: cryst_private.get_k_cells_num_parallelogram(
            arg[0], arg[0], points, v1, v2, 0, x0_equivalent, ceil(frame_scale_v1), ceil(frame_scale_v2), v1_frame, v2_frame, limits_extenders)
        gen_num_of_1_cells = lambda arg: cryst_private.get_k_cells_num_parallelogram(
            arg[0], arg[0], edges, v1, v2, 1, x0_equivalent, ceil(frame_scale_v1), ceil(frame_scale_v2), v1_frame, v2_frame, limits_extenders)
        gen_num_of_2_cells = lambda arg: cryst_private.get_k_cells_num_parallelogram(
            arg[0], arg[0], faces, v1, v2, 2, x0_equivalent,ceil(frame_scale_v1), ceil(frame_scale_v2), v1_frame, v2_frame, limits_extenders)
        # Find polynomials describing growth functions
        for s in range(N):
            args = tuple((N*i + s,) for i in range(2,5) )
            polynomials_0_cells.append(cryst_private.find_poly(gen_num_of_0_cells, args, symmetric_growth=symmetric_growth)[0])
            polynomials_1_cells.append(cryst_private.find_poly(gen_num_of_1_cells, args, symmetric_growth=symmetric_growth)[0])
            polynomials_2_cells.append(cryst_private.find_poly(gen_num_of_2_cells, args, symmetric_growth=symmetric_growth)[0])
        # Convert lists of polynomials to growth function
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
        # Find polynomials describing growth functions
        for s_k in range(K):
            for s_l in range(L):
                args = tuple( (s_k + arg[0] * K, s_l + arg[1] * L)  for arg in ((2,2), (2,3), (3,2), (3,3)))
                polynomials_0_cells.append(cryst_private.find_poly(gen_num_of_0_cells, args, symmetric_growth)[0])
                polynomials_1_cells.append(cryst_private.find_poly(gen_num_of_1_cells, args, symmetric_growth)[0])
                polynomials_2_cells.append(cryst_private.find_poly(gen_num_of_2_cells, args, symmetric_growth)[0])
        # Convert lists of polynomials to growth function
        return (gf.growth_function(polynomials_0_cells, 2, 0, (K, L)), gf.growth_function(polynomials_1_cells, 2, 1, (K, L)), gf.growth_function(polynomials_2_cells, 2, 2, (K, L)))

class Polygon:
"""
Class representing repeating motif
"""
    def __init__(self, cells):
        self.cells = list()
        self.cells.append(cells[0])
        # Vertices of cells with higher dimmensions are sorted to achive unique representation
        for i in range(1, len(cells)):
            self.cells.append(tuple(cryst_private.sort_points_in_cell(cell) for cell in cells[i]))

class Tessellation:
"""
Class representing tessellation
"""
    __FILE_PREFIXES = ("domains_parallelogram", "domains_lines_and_points")
    def __init__(self, polygon, cartesian_tessellation_vectors=None):
        self.polygon = polygon
        self.translation_vectors =tuple(tuple(row) for row in matrix.identity(len(self.polygon.cells) - 1))
        self.dim = len(self.polygon.cells) - 1
        self.cartesian_vectors = cartesian_tessellation_vectors

 
##    @param symmetric_growth
##    @param normalized Normalization of growth functions can be described as dividing functions' argument by 2.
##    @return Tuple containing topological growth functions.  Growth functions are sorted ascending by dimmension of corresponding cells.
    def get_growth_polynomials_parallelogram(self, scale_v1 = 1, scale_v2 = 1, x0 = (0,0), symmetric_growth=True):
    """Calculates crystallographical growth functions."""
        if self.dim != 2:
            raise NotImplementedError("It is not implemented yet.")
        return get_crystallographic_growth_functions(self.polygon.cells[0], self.polygon.cells[1], self.polygon.cells[2], self.translation_vectors[0],
                                                    self.translation_vectors[1], x0, frame_scale_v1 = scale_v1, frame_scale_v2 = scale_v2,
                                                   symmetric_growth=symmetric_growth)
    @param symmetric_growth
    @param normalized Normalization of growth functions can be described as dividing functions' argument by 2.
    @return Tuple containing topological growth functions.  Growth functions are sorted ascending by dimmension of corresponding cells.
    def get_growth_polynomials(self, symmetric_growth=True, normalized = True):
    """ Calculates topological growth functions. """
        return get_topological_growth_polynomials(self.polygon.cells, self.translation_vectors, symmetric_growth, normalized)
    """
    This method plots edges of the tessellation
    @param repetition_of_polygon Tuple of tuples describing how many times image should be repeated in each direction.
    @return Plot of tessellation
    """
    def plot_edges(self, repetition_of_polygon, thickness = 5):
        if repetition_of_polygon == None:
            repetition_of_polygon = ((0,0), (0,0))
        coord_matrix = matrix(self.cartesian_vectors)
        get_edge_in_cartesian_coordinates = lambda edge: (vector(edge[0]) * coord_matrix, vector(edge[1]) * coord_matrix)
        G = Graphics()
        for i in range(repetition_of_polygon[0][0],  repetition_of_polygon[0][1] + 1):
            for j in range(repetition_of_polygon[1][0],  repetition_of_polygon[1][1] + 1):
                G = G + sum(line(get_edge_in_cartesian_coordinates(cryst_private.translate_face(e, (i,j))), color="black", thickness = thickness )
                            for e in self.polygon.cells[1])
        return G

##    @brief Method for plotting orphic diagrams and describing them
##
##   @param symmetric_growth True, if growth is the same in each direction, False otherwise.
##   @param full_plot=False If False, only plot of 2-d domains is created
##   @param description If True, textual description of domains is printed.
##   @param export_format Format of exported image. Supported values: ".eps", ".pdf", ".pgf", ".png", ".ps", ".svg,"
##   @param file_name_sufix Sufix added to name of file containing orphic diagram
##   @param xmin Minimal x coordinate of final image. The argument doesn't have an effect if return_object = True.
##   @param xmax Maximal x coordinate of final image. The argument doesn't have an effect if return_object = True.
##   @param ymin Minimal y coordinate of final image. The argument doesn't have an effect if return_object = True.
##   @param ymax Maximal y coordinate of final image. The argument doesn't have an effect if return_object = True.
##   @param repetition_of_unit_cells Tuple of tuples describing how many times orphic diagram should be repeated in each direction. The argument doesn't have an effect if return_object = True.
##   @param repetition_of_polygon Tuple of tuples describing how many times orphic diagram should be repeated in each direction. The argument doesn't have an effect if return_object = True.
##   @param fit_image If true edges of tessellation outside orphic diagrams are not visible. The argument doesn't have an effect if return_object = True.
##   @param return_object If True instance of class @ref RegionsDescription is returned
##   @param colors_lists Tuple consisting of two lists describing colors of plots. If list is None, then colors of corresponding plot are chosen automatically. The argument doesn't have an effect if return_object = True.
    def plot_domains(self, symmetric_growth = True, full_plot=False, description = True, export_format=None, file_name_sufix = "",
                     xmin=None, xmax=None, ymin=None, ymax=None, repetition_of_unit_cells = None, repetition_of_polygon = None, fit_image=False,
                     return_object=False, colors_lists=(None, None)):
    """This method
1) plots orphic diagrams and describe them in textual form or
2) returns instance of class @ref RegionsDescription

    The mode of operation of the function depends on the argument return_object"""
        regions_desc = cryst_visualisation.find_regions(self, symmetric_growth, full_plot=full_plot)
        if return_object:
            return regions_desc
        plots = regions_desc.get_plots(repetition_of_unit_cells=repetition_of_unit_cells, repetition_of_polygon=repetition_of_polygon,
                                       colors_lists=colors_lists, fit_image=fit_image)
        for i in range(len(plots)):
            plots[i].show(aspect_ratio=1, axes=False, xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax)
            if export_format != None:
                plots[i].save(Tessellation.__FILE_PREFIXES[i] + file_name_sufix + "." + export_format,
                              xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, axes=False)

        if description:
            regions_desc.describe()


##@param file_path Path to file describing tessellation
##@param cartesian_vectors_included It should be False if input file doesn't contain information about translation vectors 
##@param dim Dimmension of tessellation
##@param crystallographic_coordinates It should be True if vertices are given in crystallographic coordinates
##
##@return Instance of Tessellation class based on input file
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
    if not (cartesian_vectors_included or crystallographic_coordinates):
        raise ValueError("If vertices are given in Cartesian coordinates, then translation vectors must be provided.")
    # Read input file
    input_file = open(file_path, 'r')
    lines = input_file.readlines()
    input_file.close()
    # Get numbers of cells
    cells_nums = tuple(int(s) for s in re.findall(r"\d+", lines[0]))
    cells = list()
    cartesian_vectors = None
    if cartesian_vectors_included:
        cartesian_vectors = tuple(cryst_private.string_to_point(lines[i-dim]) for i in range(dim))
    if crystallographic_coordinates:
        points = tuple(cryst_private.string_to_point(lines[i])
                  for i in range(1, cells_nums[0] + 1))
    else:
        # Convert Cartesian to crystallographical coordinates
        converter = lambda p: tuple(vector(p) * ~matrix(cartesian_vectors))
        points = tuple(converter(cryst_private.string_to_point(lines[i]))
                  for i in range(1, cells_nums[0] + 1))
    # Add vertices to list of cells
    cells.append(points)
    # Add cells with higher dimmension to list of cells
    acc_cells_number = cells_nums[0] + 1
    for i in range(1, dim + 1):
        cells.append(tuple(cryst_private.string_to_cell(lines[j], points)
                  for j in range(acc_cells_number, acc_cells_number + cells_nums[i])))
        acc_cells_number = acc_cells_number + cells_nums[i]
    return Tessellation(Polygon(cells),cartesian_vectors)
