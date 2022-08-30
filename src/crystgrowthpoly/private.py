""" Functions for tesselation manipulation"""
from sage.all import *

def translate_vector(tup, vec):
    return tuple(tup[i] + vec[i] for i in range(len(tup)))

def translate_face(edge, vec):
    return tuple(translate_vector(edge[i], vec) for i in range(len(edge)))

"""
Add cells to tesselation model

@param cells_in_tesselation Cells, which already are in tesselaion
@param cells All repeating figure's cells of some dimmension
@param vec Translation vector
@param move_operator Function translating given cell by vec
@param censor  Bool function returning True if cell should be counted, False otherwise
"""
def add_cells(cells_in_tesselation, cells, vec, move_operator):
    for c in cells:
        cells_in_tesselation.add(move_operator(c, vec))

"""
Add cells to tesselation model and count number of new cells

@param cells_in_tesselation Cells, which already are in tesselaion
@param cells All repeating figure's cells of some dimmension
@param vec Translation vector
@param move_operator Function translating given cell by vec
@param censor  Bool function returning True if cell should be counted, False otherwise
"""
def add_cells_censored(cells_in_tesselation, cells, vec, move_operator, censor):
    for c in cells:
        # Translate vector
        cell = move_operator(c, vec)
        if censor(cell):
            if not cell in cells_in_tesselation:
                cells_in_tesselation.add(cell)

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
    dim = len(vectors)
    return tuple( sum( vectors[j][i] *scalars[j]  for j in range(dim))   for i in range(dim))

"""
Return difference of two vectors
"""
def subtract_vectors(minuend, subtrahend):
    return tuple(z[0] - z[1] for z in zip(minuend, subtrahend))

"""
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
Returns function for checking if point is inside parallelogram
"""
def gen_is_point_in_parallelogram(n, m, scale_v1, scale_v2, x0):
    limiter_n = n * scale_v1
    limiter_m = m * scale_v2
    def is_point_in_parallelogram(point):
        vec  = subtract_vectors(point, x0)
        return  0 <= vec[0] and 0 <= vec[1] and vec[0] <= limiter_n and vec[1] <= limiter_m
    return  is_point_in_parallelogram


"""
Returns function for checking if cell is inside parallelogram defined by
"""
def gen_is_face_in_parallelogram(n, m, scale_v1, scale_v2, x0):
    is_point_in_parallelogram = gen_is_point_in_parallelogram(n, m, scale_v1, scale_v2, x0)
    return (lambda face: all(is_point_in_parallelogram(point) for point in face))


"""
Creates matrix and return its inverse
@param gen_val Function generating value of polynomial in given points

"""
def generate_inversed_matrix(args):
    if len(args[0]) == 1:
        return ~matrix(((arg[0]**2, arg[0], 1) for arg in args))
    else:
        return ~matrix(((arg[0]*arg[1], arg[0], arg[1], 1) for arg in args))

"""
@param gen_val Function generating value of polynomial in given points

"""
def find_poly(gen_val, args):
    value_vector = column_matrix((gen_val(arg) for arg in args))
    mat = generate_inversed_matrix(args)
    return (mat * value_vector).transpose()

"""
Function calculating number of 0 cells

@param n Number of translations in v1 direction
@param m Number of translations in v1 direction
@param points List of 0-cells
@param v1  Translation vector
@param v2  Translation vector

@return Number of 0 cells in tesselation
"""
def get_0_cells_num(n,m, points, v1, v2):
    points_in_tesselation = set()
    for i in range(n):
        for j in range(m):
            # Translate figure by v1*i + v2*j, add its points to points' set
            add_cells(points_in_tesselation, points, multiply_by_scalar_and_add([v1, v2], [i, j]), translate_vector)
    return len(points_in_tesselation)

def get_k_cells_num_parallelogram(n,m, cells, v1, v2, k, x0, scale_v1, scale_v2, frame_v1, frame_v2, additional_limits):
    cells_in_tesselation = set()
    if k == 0:
        # Count 0 cells
        censor = gen_is_point_in_parallelogram(n, m, scale_v1, scale_v2, x0)
        move_operator = translate_vector
    else:
        # Count 1 or 2 cells
        censor = gen_is_face_in_parallelogram(n, m, scale_v1, scale_v2, x0)
        move_operator = translate_face
    for i in range(-additional_limits[0][0] , (n + additional_limits[0][1]) * scale_v1):
        for j in range(-additional_limits[1][0] , (m + additional_limits[0][1]) * scale_v2):
            add_cells_censored(cells_in_tesselation, cells, multiply_by_scalar_and_add([v1, v2], [i, j]), move_operator, censor)
    return len(cells_in_tesselation)

def string_to_point(input_string):
    coord_strings = input_string.split()
    return vector((sage_eval(coord_strings[0]), sage_eval(coord_strings[1])), immutable=True)

def string_to_cell(input_string, points):
    coord_strings = input_string.split()
    return tuple(points[int(num_str) - 1] for num_str in coord_strings)

"""
@brief This function check how many additional translations need to be performed.

@param x0 Anchor point of the frame
@param points Vertices in repeating motif
"""
def get_limits_extenders(x0, points):
    max_v1 = 0
    min_v1 = 0
    max_v2 = 0
    min_v2 = 0
    for p in points:
        if p[0] > max_v1:
            max_v1 = p[0]
        elif p[0] < min_v1:
            min_v1 = p[0]
            if min_v1 in ZZ:
                min_v1 = min_v1 -1
        if p[1] > max_v2:
            max_v2 = p[1]
        elif p[1] < min_v2:
            min_v2 = p[1]
            if min_v2 in ZZ:
                min_v2 = min_v2 - 1
    return ((floor(max_v1), -floor(min_v1)), (floor(max_v2), -floor(min_v2)))
        
        