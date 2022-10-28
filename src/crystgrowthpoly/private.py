""" Functions for tessellation manipulation"""
from sage.all import *
import itertools
import functools

def translate_vector(tup, vec):
    return tuple(tup[i] + vec[i] for i in range(len(tup)))

"""
Translate cell, which dimmension is higher then 1
@param face Cell to translate
@param vec Translation vector
"""
def translate_face(face, vec):
    return tuple(translate_vector(face[i], vec) for i in range(len(face)))

"""
Add translate and add cells to tessellation model

@param cells_in_tessellation Cells, which already are in tessellaion
@param cells All repeating figure's cells of some dimmension
@param vec Translation vector
@param move_operator Function translating given cell by vec
"""
def add_cells(cells_in_tessellation, cells, vec, move_operator):
    for c in cells:
        cells_in_tessellation.add(move_operator(c, vec))

"""
Add cells to tessellation model and count number of new cells which satisfy given condition

@param cells_in_tessellation Cells, which already are in tessellation
@param cells All repeating figure's cells of some dimmension
@param vec Translation vector
@param move_operator Function translating given cell by vec
@param censor  Bool function returning True if cell should be counted, False otherwise
"""
def add_cells_censored(cells_in_tessellation, cells, vec, move_operator, censor):
    for c in cells:
        cell = move_operator(c, vec)
        if censor(cell):
            cells_in_tessellation.add(cell)

"""
Returns vector multiplied by scalar
"""
def multiply_vector(vec, scalar):
    return tuple(a * scalar for a in vec)

"""
Multiplies vectors by scalars and add them.

@param vectors List of vectors
@param scalars List of scalars
"""
def multiply_by_scalar_and_add(vectors, scalars):
    dim = len(vectors[0])
    return tuple(sum(vectors[j][i] * scalars[j] for j in range(len(vectors))) for i in range(dim))

"""
Return difference of two vectors
"""
def subtract_vectors(minuend, subtrahend):
    return tuple(z[0] - z[1] for z in zip(minuend, subtrahend))

"""
Sorts points representing cell

@param cell Tuple representing cell
@returns Sorted tuple representing input cell
"""
def sort_points_in_cell(cell):
    def is_grater_then(tup1, tup2):
        for i in range(len(tup1)):
            if tup1[i] > tup2[i]:
                return true
        return false
    cell_to_list = list(cell)
    cell_to_list.sort()
    return tuple(cell_to_list)

"""
Returns function checking if point is inside parallelogram

@param n Length of rectangle's side before scaling (In crystallographical coordinates parallelogram frames used by us are rectangles)
@param m Lenght of rectangle's side before scaling (In crystallographical coordinates parallelogram frames used by us are rectangles)
@param scale_v1 Scaling value
@param scale_v2 Scaling value
@param x0 Anchor point of parallelogram frame
"""
def gen_is_point_in_parallelogram(n, m, scale_v1, scale_v2, x0):
    limiter_n = n * scale_v1
    limiter_m = m * scale_v2
    def is_point_in_parallelogram(point):
        vec  = subtract_vectors(point, x0)
        return  0 <= vec[0] and 0 <= vec[1] and vec[0] <= limiter_n and vec[1] <= limiter_m
    return  is_point_in_parallelogram


"""
Returns function checking if cell (with dimmension higher then 1) is inside parallelogram

@param n Length of rectangle's side before scaling (In crystallographical coordinates parallelogram frames used by us are rectangles)
@param m Lenght of rectangle's side before scaling (In crystallographical coordinates parallelogram frames used by us are rectangles)
@param scale_v1 Scaling value
@param scale_v2 Scaling value
@param x0 Anchor point of parallelogram frame
"""
def gen_is_face_in_parallelogram(n, m, scale_v1, scale_v2, x0):
    is_point_in_parallelogram = gen_is_point_in_parallelogram(n, m, scale_v1, scale_v2, x0)
    return (lambda face: all(is_point_in_parallelogram(point) for point in face))


"""
Creates matrix used for growth polynomial interpolation

@param args Data points for interpolation
@param symmetric_growth True, if growth is the same in each direction, False otherwise.
"""
def generate_inversed_matrix(args, symmetric_growth):
    if symmetric_growth:
        n = len(args) -1
        return ~matrix(tuple(tuple(arg[0]**(n - i) for i in range(len(args))) for arg in args))
    else:
        dim = len(args[0])
        indexes  = range(dim)
        row_scheme = tuple(itertools.chain(*(tuple(itertools.combinations(indexes, dim - n)) for n in range(dim + 1))))
        def generate_row(arg):
            return tuple(functools.reduce(lambda acc, index,:  acc * arg[index], chosen_indexes, 1) for chosen_indexes in row_scheme)
        return ~matrix(generate_row(arg) for arg in args)

"""
@param gen_val Function generating value of polynomial in given points
@param args Data points
@param symmetric_growth True, if growth is the same in each direction, False otherwise.

@return Coefficients of polynomial
"""
def find_poly(gen_val, args, symmetric_growth):
    value_vector = column_matrix(tuple(gen_val(arg) for arg in args))
    mat = generate_inversed_matrix(args, symmetric_growth)
    return (mat * value_vector).transpose()

"""
Function calculating number of cells with some dimmension after arguments[i] steps in direction translation_vectors[i]

@param arguments Numbers of steps
@param cells List of cells with given dimmension in repeating unit
@param translation_vectors List of translation vectors
@return Number of cells with some dimmension after arguments[i] steps in direction translation_vectors[i]
"""
def get_k_cells_num(arguments, cells, translation_vectors):
    cells_in_tessellation = set()
    if type(cells[0][0]) != tuple:
        # Count 0 cells
        move_operator = translate_vector
    else:
        # Count higher dimmension cells
        move_operator = translate_face
    for translations_numbers in itertools.product(*(range(arg) for arg in arguments)):
        add_cells(cells_in_tessellation, cells, multiply_by_scalar_and_add(translation_vectors, translations_numbers), move_operator)
    return len(cells_in_tessellation)

"""
Function calculating number of cells with given dimmension inside parallelogram

@param n Number of translations in v1 direction
@param m Number of translations in v2 direction
@param cells List of k-cells in repeating unit
@param v1  Translation vector
@param v2  Translation vector
@param k Dimmension of cell
@param x0 Anchor point of parallelogram frame
@param scale_v1 Scaling value for vector v1
@param scale_v2 Scaling value for vector v2
@param frame_v1 v1 * scale_v1
@param frame_v2 v2 * scale_v2
@param additional_limits numbers of additional translations that need to be performed

@return Number of k-cells inside parallelogram
"""
def get_k_cells_num_parallelogram(n,m, cells, v1, v2, k, x0, scale_v1, scale_v2, frame_v1, frame_v2, additional_limits):
    cells_in_tessellation = set()
    if k == 0:
        # Count 0 cells
        censor = gen_is_point_in_parallelogram(n, m, scale_v1, scale_v2, x0)
        move_operator = translate_vector
    else:
        # Count higher dimmension cells
        censor = gen_is_face_in_parallelogram(n, m, scale_v1, scale_v2, x0)
        move_operator = translate_face
    for i in range(-additional_limits[0][0] , (1 + n + additional_limits[0][1]) * scale_v1):
        for j in range(-additional_limits[1][0] , (1 + m + additional_limits[1][1]) * scale_v2):
            add_cells_censored(cells_in_tessellation, cells, multiply_by_scalar_and_add([v1, v2], [i, j]), move_operator, censor)
    return len(cells_in_tessellation)

def string_to_point(input_string):
    coord_strings = input_string.split()
    return tuple(sage_eval(string) for string in coord_strings)

"""
@brief Converts string to cell

@param x0 input_string String encoding cell consisting of vertices' indexes
@param points Vertices in repeating motif
"""
def string_to_cell(input_string, points):
    coord_strings = input_string.split()
    return tuple(points[int(num_str) - 1] for num_str in coord_strings)

"""
@brief This function calculate how many additional translations are needed to be performed.

@param x0 Anchor point of the frame
@param points Vertices in repeating motif
"""
def get_limits_extenders(x0, points):
    max_v1 = 0
    min_v1 = 0
    max_v2 = 0
    min_v2 = 0
    # Find extreme coordinates of repeating unit
    for p in points:
        if p[0] - x0[0] > max_v1:
            max_v1 = p[0] - x0[0]
        elif p[0] - x0[0] < min_v1:
            min_v1 = p[0] - x0[0]
            if min_v1 in ZZ:
                min_v1 = min_v1 -1
        if p[1] - x0[1] > max_v2:
            max_v2 = p[1] - x0[1]
        elif p[1] - x0[1] < min_v2:
            min_v2 = p[1] - x0[1]
            if min_v2 in ZZ:
                min_v2 = min_v2 - 1
    return ((floor(max_v1), -floor(min_v1)), (floor(max_v2), -floor(min_v2)))

"""
This function calculates alternative sums of corresponding coefficients.

@param coefficient_lists List of lists containing polynomials' coefficients
@return List of alternative sums
"""
def get_alternative_sum(coefficient_lists):
    length_of_lists = len(coefficient_lists[0]) # Number of single polynomial's coefficients
    return_list = list()
    for i in range(length_of_lists):
        acc = 0
        for j  in range(len(coefficient_lists)):
            acc = acc + coefficient_lists[j][i] * ((-1) ** (j % 2))
        return_list.append(acc)
    return return_list

"""
Calculates coefficients of polynomial corresponding to a cell with dimension dim - 1

@param coefficient_lists List containing lists of growth polynomial coefficients corresponding to a cells with dimensions 0 through dim - 2 
@param max_dimmension_coefficients Coefficients of growth polynomial corresponding to cells of dimmension dim
@param dim Dimmension of tessellation
@return Coefficients of growth polynomial corresponding to cells with dimmension dim - 1
"""
def calculate_coefficients_based_ne_euler_charcteristics(coefficient_lists, max_dimmension_coefficients, dim):
    return_list = get_alternative_sum(coefficient_lists)
    return_list[-1] = return_list[-1] - 1
    if dim % 2 == 0:
        for i in range(len(return_list)):
            return_list[i] = return_list[i] + max_dimmension_coefficients[i]
    else:
        for i in range(len(return_list)):
            return_list[i] = -(return_list[i] - max_dimmension_coefficients[i])
    return return_list