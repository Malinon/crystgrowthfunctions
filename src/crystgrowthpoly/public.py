import crystgrowthpoly.private as cryst_private
import crystgrowthpoly.visualisation as cryst_visualisation
import crystgrowthpoly.growth_function as gf
import re
from sage.all import *

"""
Function finding growth polynomials

@param points 0-cells of repeating figure
@param num_of_2_cells Number of 2-cells in repeating figure
@param v1  tesselation vector
@param v2  tesselation vector
"""
def get_topological_growth_polynomials(points, num_of_2_cells, v1, v2, symmetric_frame=True):
    # Generate function calculating number of 0-cells in given step of tesselation
    if symmetric_frame:
        args = (2, 3, 4)
        variables_num = 1
        gen_val = lambda arg: cryst_private.get_0_cells_num(arg, arg, points, v1, v2)
        polynomial_0_cells_coeff = cryst_private.find_poly(gen_val, args)[0]
        polynomial_2_cells_coeff = (num_of_2_cells, 0, 0)
        polynomial_1_cells_coeff = (polynomial_0_cells_coeff[0] + num_of_2_cells, polynomial_0_cells_coeff[1], polynomial_0_cells_coeff[2] - 1)
    else:
        gen_val = lambda arg: cryst_private.get_0_cells_num(arg[0], arg[1], points, v1, v2)
        args = ((2,2), (2,3), (3,2), (3,3))
        variables_num = 2
        polynomial_0_cells_coeff = cryst_private.find_poly_2_variables(gen_val, args)[0]
        polynomial_2_cells_coeff = (num_of_2_cells, 0, 0,0)
        polynomial_1_cells_coeff = (polynomial_0_cells_coeff[0] + num_of_2_cells, polynomial_0_cells_coeff[1],
                                    polynomial_0_cells_coeff[2], polynomial_0_cells_coeff[3] - 1)
    return (gf.growth_function((polynomial_0_cells_coeff,), variables_num, 0),
            gf.growth_function((polynomial_1_cells_coeff,), variables_num, 1),
            gf.growth_function((polynomial_2_cells_coeff,), variables_num, 2))


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
    polynomials_0_cells = []
    polynomials_1_cells = []
    polynomials_2_cells = []
    if in_one_variable:
        N = lcm(rational_scale_v1.denominator(), rational_scale_v2.denominator())
        if N == 1:
            polynomial_finder = cryst_private.find_poly_par
        else:
            polynomial_finder = cryst_private.find_poly
        # Prepare functions calculating numbers of cells
        gen_num_of_0_cells = lambda arg: cryst_private.get_k_cells_num_parallelogram(
            arg[0], arg[0], points, v1, v2, 0, x0_equivalent, ceil(frame_scale_v1), ceil(frame_scale_v2), v1_frame, v2_frame)
        gen_num_of_1_cells = lambda arg: cryst_private.get_k_cells_num_parallelogram(
            arg[0], arg[0], edges, v1, v2, 1, x0_equivalent, ceil(frame_scale_v1), ceil(frame_scale_v2), v1_frame, v2_frame)
        gen_num_of_2_cells = lambda arg: cryst_private.get_k_cells_num_parallelogram(
            arg[0], arg[0], faces, v1, v2, 2, x0_equivalent,ceil(frame_scale_v1), ceil(frame_scale_v2), v1_frame, v2_frame)
        for s in range(N):
            args = tuple((N*i + s,) for i in range(2,5) )
            # Find growth polynomials
            polynomials_0_cells.append(polynomial_finder(gen_num_of_0_cells, args)[0])
            polynomials_1_cells.append(polynomial_finder(gen_num_of_1_cells, args)[0])
            polynomials_2_cells.append(polynomial_finder(gen_num_of_2_cells, args)[0])
        return (gf.growth_function(polynomials_0_cells, 1, 0, N), gf.growth_function(polynomials_1_cells, 1, 1, N), gf.growth_function(polynomials_2_cells, 1, 2, N))
    else:
        K = rational_scale_v1.denominator()
        L = rational_scale_v2.denominator()
        if K == 1 and L == 1:
            polynomial_finder = cryst_private.find_poly_par
        else:
            polynomial_finder = cryst_private.find_poly_2_variables
        gen_num_of_0_cells = lambda arg: cryst_private.get_k_cells_num_parallelogram(
            arg[0], arg[1], points, v1, v2, 0, x0_equivalent, ceil(frame_scale_v1), ceil(frame_scale_v2), v1_frame, v2_frame)
        gen_num_of_1_cells = lambda arg: cryst_private.get_k_cells_num_parallelogram(
            arg[0], arg[1], edges, v1, v2, 1, x0_equivalent, ceil(frame_scale_v1), ceil(frame_scale_v2), v1_frame, v2_frame)
        gen_num_of_2_cells = lambda arg: cryst_private.get_k_cells_num_parallelogram(
            arg[0], arg[1], faces, v1, v2, 2, x0_equivalent,ceil(frame_scale_v1), ceil(frame_scale_v2), v1_frame, v2_frame)
        for s_k in range(K):
            for s_l in range(L):
                args = tuple( (s_k + arg[0] * K, s_l + arg[1] * L)  for arg in ((2,2), (2,3), (3,2), (3,3)))
                polynomials_0_cells.append(polynomial_finder(gen_num_of_0_cells, args)[0])
                polynomials_1_cells.append(polynomial_finder(gen_num_of_1_cells, args)[0])
                polynomials_2_cells.append(polynomial_finder(gen_num_of_2_cells, args)[0])
        return (gf.growth_function(polynomials_0_cells, 2, 0, (K, L)), gf.growth_function(polynomials_1_cells, 2, 1, (K, L)), gf.growth_function(polynomials_2_cells, 2, 2, (K, L)))

class Polygon:
    def __init__(self,points, edges, faces):
        self.points = points
        self.edges = tuple(cryst_private.sort_points_in_cell(edge) for edge in edges)
        self.faces = tuple(cryst_private.sort_points_in_cell(face) for face in faces)

class Tesselation:
    def __init__(self, polygon, v1, v2, cartesian_tesselation_vectors=None):
        self.polygon = polygon
        self.v1 = v1
        self.v2 = v2
        self.cartesian_vectors = cartesian_tesselation_vectors
    def get_growth_polynomials_parallelogram(self, scale_v1 = 1, scale_v2 = 1, x0 = (0,0), symmetric_frame=True):
        return get_crystalographic_growth_functions(self.polygon.points, self.polygon.edges, self.polygon.faces, self.v1, self.v2, x0,
                                                    frame_scale_v1 = scale_v1, frame_scale_v2 = scale_v2,
                                                   in_one_variable=symmetric_frame)
    def get_growth_polynomials(self, symmetric_frame=True):
        return get_topological_growth_polynomials(self.polygon.points, len(self.polygon.faces), self.v1, self.v2, symmetric_frame)
    def plot_edges(self):
        coord_matrix = matrix(self.cartesian_vectors)
        get_edge_in_cartesian_coordinates = lambda edge: (vector(edge[0]) * coord_matrix, vector(edge[1]) * coord_matrix)
        return sum( line(get_edge_in_cartesian_coordinates(e), color="black") for e in self.polygon.edges)
    def plot_domains(self, symmetric_growth = True, full_plot=False, description = True):
        regions = cryst_visualisation.find_regions(self, symmetric_growth, full_plot=full_plot)
        own_plot = self.plot_edges()
        (regions.parallelograms_plot + own_plot).show(aspect_ratio=1, axes=False)
        if full_plot:
            (own_plot + regions.edges_plot + regions.corners_plot).show(aspect_ratio=1, axes=False)
        if description:
            for polynomials in regions.regions_description.keys():
                print("***************************")
                print("For x_0 in")
                for domain in regions.regions_description[polynomials]:
                    domain.describe(self.v1, self.v2)
                for f in regions.regions_description[polynomials][0].growth_f:
                    f.show()
                print("***************************")


"""
If crystalographic_coordinates is set True, then cartesian_vectors_included is ignored
"""
def read_tessellation_from_file(file_path, crystalographic_coordinates=True, cartesian_vectors_included=True):
    input_file = open(file_path, 'r')
    lines = input_file.readlines()
    cells_nums = tuple(int(s) for s in re.findall(r"\d+", lines[0]))
    cartesian_v1 = None
    cartesian_v2 = None
    if cartesian_vectors_included:
        cartesian_v1 = cryst_private.string_to_point(lines[ cells_nums[0] + cells_nums[1] + cells_nums[2] + 1])
        cartesian_v2 = cryst_private.string_to_point(lines[ cells_nums[0] + cells_nums[1] + cells_nums[2] + 2])
    if crystalographic_coordinates:
        points = [cryst_private.string_to_point(lines[i]) for i in range(1, cells_nums[0] + 1)]
        v1 = (1, 0)
        v2 = (0, 1)
    else:
        v1 = cartesian_v1
        v2 = cartesian_v2
        points = tuple(multiply_by_scalar_and_add([v1, v2], string_to_point(lines[i])) for i in range(1, cells_nums[0] + 1))
        cartesian_vectors = (v1, v2)
    edges = tuple(cryst_private.string_to_cell(lines[i], points) for i in range(1+cells_nums[0], cells_nums[0] + cells_nums[1] + 1))
    faces = tuple(cryst_private.string_to_cell(lines[i], points) for i in range(1+cells_nums[0] + cells_nums[1], cells_nums[0] + cells_nums[1] + cells_nums[2] + 1))
    return Tesselation(Polygon(points, edges, faces), v1, v2, (cartesian_v1, cartesian_v2))
