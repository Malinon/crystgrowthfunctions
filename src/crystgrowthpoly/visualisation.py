import crystgrowthpoly.private as cryst_private
from sage.all import *
import itertools

def find_mantisse(tes, coordinate, is_vertical):
    if is_vertical:
        index = 0
    else:
        index = 1
    min_mantisse = 1
    for p in tes:
        m = calculate_mantisse(p[index] -  coordinate)
        if m < min_mantisse and m != 0:
            min_mantisse = m

    return min(coordinate + min_mantisse, 1)


def find_boundary_lines(tes, is_vertical):
    coordinate = 0
    endpoints = [coordinate]
    while  coordinate != 1:
        coordinate = find_mantisse(tes, coordinate, is_vertical)
        endpoints.append(coordinate)
    return endpoints

class OpenRectangle:
    def __init__(self, lower_left_corner, higher_right_corner):
        self.higher_right_corner = higher_right_corner
        self.lower_left_corner = lower_left_corner
        self.cartesian_corners = None
    def get_point_from_interior(self):
        return ((self.higher_right_corner[0]+self.lower_left_corner[0])/Integer(2), (self.lower_left_corner[1]+self.higher_right_corner[1])/Integer(2) )
    def set_growth_function(self, tes, symmetric_frame):
        self.growth_f = tes.get_growth_polynomials_parallelogram(1,1,self.get_point_from_interior(), symmetric_frame)
    def get_polynomials():
        return(growth_f[0].polynomials, growth_f[1].polynomials, growth_f[2].polynomials)
    def describe(self, v1, v2):
        p2 = vector((self.higher_right_corner[0], self.lower_left_corner[1]))
        p4 = vector((self.lower_left_corner[0], self.higher_right_corner[1]))
        print("Parallelogram ", (self.lower_left_corner, p2, self.higher_right_corner, p4))
    def to_graphic_polygon(self, v1, v2, color, translation=(0,0)):
        if self.cartesian_corners == None:
            p2 = vector((self.higher_right_corner[0], self.lower_left_corner[1]))
            p4 = vector((self.lower_left_corner[0], self.higher_right_corner[1]))
            inv_mat = column_matrix((v1,v2))
            real_p1 = inv_mat * vector(self.lower_left_corner)
            real_p2 = inv_mat * p2
            real_p3 = inv_mat * vector(self.higher_right_corner)
            real_p4 = inv_mat * p4
            self.cartesian_corners = (real_p1, real_p2, real_p3, real_p4)
        return polygon((corner + vector(translation)  for corner in self.cartesian_corners), color=color)
    def get_polynomials(self):
        return(self.growth_f[0].polynomials, self.growth_f[1].polynomials, self.growth_f[2].polynomials)

class OpenLine:
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.cartesian_endpoints = None
    def plot(self, v1, v2, color, translation=(0,0)):
        if self.cartesian_endpoints == None:
            mat = matrix((v1, v2))
            real_start = vector(self.start)  * mat
            real_end = vector(self.end)  * mat
            self.cartesian_endpoints = (real_start, real_end)
        return line((self.cartesian_endpoints[0] + vector(translation), self.cartesian_endpoints[1] + vector(translation)), color=color)
    def get_point_from_interior(self):
        return ((self.start[0] + self.end[0]) / Integer(2), (self.start[1] + self.end[1]) / Integer(2))
    def get_polynomials():
        return(growth_f[0].polynomials, growth_f[1].polynomials, growth_f[2].polynomials)
    def set_growth_function(self, tes, symmetric_frame):
        self.growth_f = tes.get_growth_polynomials_parallelogram(1,1,self.get_point_from_interior(), symmetric_frame)
    def describe(self, v1, v2):
        print("Line ", (self.start, self.end))
    def get_polynomials(self):
        return(self.growth_f[0].polynomials, self.growth_f[1].polynomials, self.growth_f[2].polynomials)

class SpecialPoint:
    def __init__(self, point):
        self.point = point
        self.cartesian_point = None
    def plot(self, v1, v2, color, translation):
        if self.cartesian_point == None:
            self.cartesian_point = vector(self.point) * matrix((v1, v2))
        return point(self.cartesian_point + vector(translation), color =color, size = 20, zorder=2)
    def get_point(self):
        return self.point
    def set_growth_function(self, tes, symmetric_frame):
        self.growth_f = tes.get_growth_polynomials_parallelogram(1,1,self.get_point(), symmetric_frame)
    def describe(self, v1, v2):
        print("Point ", self.point)
    def get_polynomials(self):
        return(self.growth_f[0].polynomials, self.growth_f[1].polynomials, self.growth_f[2].polynomials)

class RegionsDescription:
    def __init__(self, parallelograms_plot, edges_plot, corners_plot, regions_description):
        self.parallelograms_plot = parallelograms_plot
        self.edges_plot = edges_plot
        self.corners_plot = corners_plot
        self.regions_description = regions_description

def calculate_mantisse(x):
    return x - floor(x)

# I only need to consider changes in points. Cell changes frames iff one f its points changes frame


def find_regions(tessellation, symmetric_frame=True, full_plot=False, repetition_of_unit_cells=None):
    points = tessellation.polygon.points
    tes = set()
    # Create fragment of tessellation
    for k in range(-1, 4):
        for l in range(-1, 4):
            cryst_private.add_cells(tes, points, cryst_private.multiply_by_scalar_and_add([(1,0), (0,1)], [k, l]), cryst_private.translate_vector)
    x_endpoints = find_boundary_lines(tes, True)
    y_endpoints = find_boundary_lines(tes, False)

    found_rectangles = tuple( OpenRectangle((x_endpoints[x_iter - 1], y_endpoints[y_iter - 1]), (x_endpoints[x_iter], y_endpoints[y_iter]))
                    for x_iter in range(1, len(x_endpoints) )
                    for y_iter in range(1, len(y_endpoints) ))
    #print("We have found rectangles!!!")
    points_list = tuple(SpecialPoint((x, y)) for x in x_endpoints for y in y_endpoints)
    #print("We have found special points!!!")
    vertical_edges_iterator = (OpenLine( (x, y_endpoints[y_iter - 1]), (x, y_endpoints[y_iter]))  for y_iter in range(1, len(y_endpoints)) for x in x_endpoints)
    horisontal_edges_iterator = (OpenLine( (x_endpoints[x_iter - 1], y), (x_endpoints[x_iter], y))  for x_iter in range(1, len(x_endpoints)) for y in y_endpoints)
    lines_list = tuple(l for l in itertools.chain(vertical_edges_iterator, horisontal_edges_iterator))
    different_types = set()
    #print("We have found lists!!!")
    polynomials = set()
    for rec in found_rectangles:
        rec.set_growth_function(tessellation, symmetric_frame)
        polynomials.add(rec.get_polynomials())
    print("Number of parallelograms colors", len(polynomials))
    colors = rainbow(len(polynomials))
    d = dict(zip(polynomials, colors))
        #print(rec.get_random_point()).translate
    if repetition_of_unit_cells == None:
        max_x = max(tessellation.polygon.points, key = lambda p: p[0])[0]
        min_x = min(tessellation.polygon.points, key = lambda p: p[0])[0]
        min_y = min(tessellation.polygon.points, key = lambda p: p[1])[1]
        max_y = max(tessellation.polygon.points, key = lambda p: p[1])[1]
    else:
        min_x = repetition_of_unit_cells[0][0]
        max_x = repetition_of_unit_cells[0][1]
        min_y = repetition_of_unit_cells[1][0]
        max_y = repetition_of_unit_cells[1][1]
    G = Graphics()
    for k in range(floor(min_x), ceil(max_x)):
        for l in range(floor(min_y), ceil(max_y)):
            G = G +  sum(rect.to_graphic_polygon(tessellation.cartesian_vectors[0], tessellation.cartesian_vectors[1], d[rect.get_polynomials()],
                                                                cryst_private.multiply_by_scalar_and_add(tessellation.cartesian_vectors, (k,l)) ) for rect in found_rectangles)
    #print("Main plot is ready")
    if full_plot:
        G_lines = Graphics()
        G_points = Graphics()
        polynomials_lines_points = set()
        for l in lines_list:
            l.set_growth_function(tessellation, symmetric_frame)
            polynomials_lines_points.add(l.get_polynomials())
        #print("We have functions for edges")
        for p in points_list:
            p.set_growth_function(tessellation, symmetric_frame)
            polynomials_lines_points.add(p.get_polynomials())
        #print("We have functions for points")
        d2 = dict(zip(polynomials_lines_points, rainbow(len(polynomials_lines_points))))
        for k in range(floor(min_x), ceil(max_x)):
            for l in range(floor(min_y), ceil(max_y)):
                G_lines = G_lines +  sum(line.plot(tessellation.cartesian_vectors[0], tessellation.cartesian_vectors[1], d2[line.get_polynomials()],
                                                                cryst_private.multiply_by_scalar_and_add(tessellation.cartesian_vectors, (k,l)) ) for line in lines_list)
                G_points = G_points +  sum(p.plot(tessellation.cartesian_vectors[0], tessellation.cartesian_vectors[1], d2[p.get_polynomials()],
                                                                cryst_private.multiply_by_scalar_and_add(tessellation.cartesian_vectors, (k,l)) ) for p in points_list)
        for poly in polynomials_lines_points:
            polynomials.add(poly)
    description_dic = dict( (polys, []) for polys in polynomials )
    for rect in found_rectangles:
        description_dic[rect.get_polynomials()].append(rect)
    if full_plot:
        for l in lines_list:
            description_dic[l.get_polynomials()].append(l)
        for p in points_list:
            description_dic[p.get_polynomials()].append(p)
        return RegionsDescription(G, G_lines, G_points,  description_dic)
    return RegionsDescription(G, None, None,  description_dic)