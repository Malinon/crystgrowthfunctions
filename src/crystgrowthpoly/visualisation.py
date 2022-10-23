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
    def __hash__(self):
        return hash((self.start, self.end))
    def __eq__(self, obj):
        return self.start == obj.start and self.end == obj.end
    def plot(self, v1, v2, color, translation=(0,0), thickness=5):
        if self.cartesian_endpoints == None:
            mat = matrix((v1, v2))
            real_start = vector(self.start)  * mat
            real_end = vector(self.end)  * mat
            self.cartesian_endpoints = (real_start, real_end)
        return line((self.cartesian_endpoints[0] + vector(translation), self.cartesian_endpoints[1] + vector(translation)), color=color, thickness=thickness)
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
    def plot(self, v1, v2, color, translation, size):
        if self.cartesian_point == None:
            self.cartesian_point = vector(self.point) * matrix((v1, v2))
        return point(self.cartesian_point + vector(translation), color =color, size = size, zorder=2)
    def get_point(self):
        return self.point
    def set_growth_function(self, tes, symmetric_frame):
        self.growth_f = tes.get_growth_polynomials_parallelogram(1,1,self.get_point(), symmetric_frame)
    def describe(self, v1, v2):
        print("Point ", self.point)
    def get_polynomials(self):
        return(self.growth_f[0].polynomials, self.growth_f[1].polynomials, self.growth_f[2].polynomials)

class RegionsDescription:
    def __init__(self, regions_lists, dimmension, tessellation):
        self.regions_lists = regions_lists
        self.__dimmension = dimmension
        self.tessellation = tessellation
        # Group regions by polynomials
        polynomials_all = set(reg.get_polynomials()
                              for reg_list in self.regions_lists
                              for reg in reg_list)
        self.description_dict = dict( (poly, []) for poly in polynomials_all)
        for reg_list in self.regions_lists:
            for reg in reg_list:
                self.description_dict[reg.get_polynomials()].append(reg)
        if dimmension == 2:
            # We use dictionaries to preserve order
            self.polynomials_par = dict((reg.get_polynomials(), None) for reg in self.regions_lists[0])
            if len(self.regions_lists) > 1:
                self.polynomials_lines_points = dict((reg.get_polynomials(), None)
                                                   for ind in range(1, len(self.regions_lists)) for reg in self.regions_lists[ind])
    def get_colors_num(self):
        if len(self.regions_lists) > 1:
            return (len(self.polynomials_par), len(self.polynomials_par))
        else:
            return len(self.polynomials_lines_points)
    def generate_censoring_polygon(self, plot_basic, repetition_of_unit_cells):
        polygon_vertice1 = cryst_private.multiply_by_scalar_and_add(self.tessellation.cartesian_vectors,
                                                               (repetition_of_unit_cells[0][0], repetition_of_unit_cells[1][0]))
        polygon_vertice2 = cryst_private.translate_vector(polygon_vertice1, cryst_private.multiply_vector(self.tessellation.cartesian_vectors[1],
                                                                                                    repetition_of_unit_cells[1][1] - repetition_of_unit_cells[1][0]))
        margin = max(abs(plot_basic.xmin() - plot_basic.xmax()), abs(plot_basic.ymin() - plot_basic.ymax())) * 0.01
        censoring_polygon = polygon(((plot_basic.xmin() - margin, plot_basic.ymin() - margin),
                            (plot_basic.xmax()  +margin, plot_basic.ymin() - margin),
                            (plot_basic.xmax()  + margin, plot_basic.ymax()  + margin),
                            (plot_basic.xmin() - margin, plot_basic.ymax() + margin),
                            polygon_vertice2,
                            cryst_private.multiply_by_scalar_and_add((polygon_vertice1, self.tessellation.cartesian_vectors[0], self.tessellation.cartesian_vectors[1]),
                                                                    (1, repetition_of_unit_cells[0][1] - repetition_of_unit_cells[0][0],
                                                                     repetition_of_unit_cells[1][1] - repetition_of_unit_cells[1][0]) ),
                            cryst_private.translate_vector(polygon_vertice1, cryst_private.multiply_vector(self.tessellation.cartesian_vectors[0],
                                                                                                    repetition_of_unit_cells[0][1] - repetition_of_unit_cells[0][0])),
                            polygon_vertice1,
                            polygon_vertice2, (plot_basic.xmin() - margin, plot_basic.ymax() + margin)), color="white", zorder=5)
        return censoring_polygon

    """
    @param repetition_of_unit_cells Pair of pairs defining number of copies of unit cells in each directions
    @param colors_lists pair containing list of colors (First for 2-d regions, second for 1-d, 0-d regions)
    @param repetition_of_polygon Pair of pairs defining number of copies of polygon in each directions
    @param fit_image If true, then polygon drawings outside colored unit cells are not visible
    """
    def get_plots(self, repetition_of_unit_cells=None, colors_lists=(None, None),
                  repetition_of_polygon=None, fit_image=False, lines_thickness = 5, point_size = 20):
        if self.__dimmension != 2:
            # It is only for 2-D
            pass
        else:
            if repetition_of_unit_cells == None:
                max_x = max(self.tessellation.polygon.cells[0], key = lambda p: p[0])[0]
                min_x = min(self.tessellation.polygon.cells[0], key = lambda p: p[0])[0]
                min_y = min(self.tessellation.polygon.cells[0], key = lambda p: p[1])[1]
                max_y = max(self.tessellation.polygon.cells[0], key = lambda p: p[1])[1]
            else:
                min_x = repetition_of_unit_cells[0][0]
                max_x = repetition_of_unit_cells[0][1]
                min_y = repetition_of_unit_cells[1][0]
                max_y = repetition_of_unit_cells[1][1]
            if colors_lists[0] == None:
                # If user hasn't selected colors, then do it automaticly
                colors_poly_dict = dict(zip(self.polynomials_par.keys(), rainbow(len(self.polynomials_par))))
            else:
                colors_poly_dict = dict(zip(self.polynomials_par.keys(), colors_lists[0]))
            G = Graphics()
            for k in range(floor(min_x), ceil(max_x)):
                for l in range(floor(min_y), ceil(max_y)):
                    G = G +  sum(rect.to_graphic_polygon(
                            self.tessellation.cartesian_vectors[0], self.tessellation.cartesian_vectors[1], colors_poly_dict[rect.get_polynomials()],
                            cryst_private.multiply_by_scalar_and_add(self.tessellation.cartesian_vectors, (k,l)) ) for rect in self.regions_lists[0])
            tes_plot = self.tessellation.plot_edges(repetition_of_polygon)
            G = G + tes_plot
            if fit_image:
                censoring_polygon = self.generate_censoring_polygon(G, repetition_of_unit_cells)
                G = G + censoring_polygon
            if len(self.regions_lists) > 1:
                G_lines = Graphics()
                G_points = Graphics()
                if colors_lists[0] == None:
                    # If user hasn't selected colors, then do it automaticly
                    colors_poly_lines_points_dict = dict(zip(self.polynomials_lines_points.keys(), rainbow(len(self.polynomials_lines_points))))
                else:
                    colors_poly_lines_points_dict = dict(zip(self.polynomials_lines_points.keys(), colors_lists[1]))
                for k in range(floor(min_x), ceil(max_x)):
                    for l in range(floor(min_y), ceil(max_y)):
                        G_lines = G_lines + sum(rect.plot(
                            self.tessellation.cartesian_vectors[0], self.tessellation.cartesian_vectors[1], colors_poly_lines_points_dict[rect.get_polynomials()],
                            cryst_private.multiply_by_scalar_and_add(self.tessellation.cartesian_vectors, (k,l)), thickness = lines_thickness) for rect in self.regions_lists[1])
                        G_points = G_points + sum(rect.plot(
                            self.tessellation.cartesian_vectors[0], self.tessellation.cartesian_vectors[1], colors_poly_lines_points_dict[rect.get_polynomials()],
                            cryst_private.multiply_by_scalar_and_add(self.tessellation.cartesian_vectors, (k,l)), size = point_size) for rect in self.regions_lists[2])
                return (G, G_lines + G_points + tes_plot)
            else:
                return (G,)
    def describe(self):
        for polynomials in self.description_dict.keys():
            print("***************************")
            print("For x_0 in")
            for domain in self.description_dict[polynomials]:
                domain.describe(self.tessellation.v1, self.tessellation.v2)
            for f in self.description_dict[polynomials][0].growth_f:
                f.show()
            print("***************************")

def calculate_mantisse(x):
    return x - floor(x)


def find_regions(tessellation, symmetric_frame=True, full_plot=False):
    points = tessellation.polygon.cells[0]
    regions_lists = [] # Based on this list RegionsDescription object is constructed
    tes = set()
    # Create fragment of tessellation
    for k in range(-1, 4):
        for l in range(-1, 4):
            cryst_private.add_cells(tes, points, cryst_private.multiply_by_scalar_and_add([(1,0), (0,1)], [k, l]), cryst_private.translate_vector)
    x_endpoints = find_boundary_lines(tes, True)
    y_endpoints = find_boundary_lines(tes, False)

    # Find parallelogram regions
    found_rectangles = tuple(OpenRectangle((x_endpoints[x_iter - 1], y_endpoints[y_iter - 1]), (x_endpoints[x_iter], y_endpoints[y_iter]))
                             for x_iter in range(1, len(x_endpoints))
                             for y_iter in range(1, len(y_endpoints)))
    for rec in found_rectangles:
        rec.set_growth_function(tessellation, symmetric_frame)
    regions_lists.append(found_rectangles)
    # If requested find 0,1-dimmensional regions
    if full_plot:
        points_list = tuple(SpecialPoint((x, y)) for x in x_endpoints for y in y_endpoints)
        vertical_edges_iterator = (OpenLine( (x, y_endpoints[y_iter - 1]), (x, y_endpoints[y_iter]))  for y_iter in range(1, len(y_endpoints)) for x in x_endpoints)
        horisontal_edges_iterator = (OpenLine( (x_endpoints[x_iter - 1], y), (x_endpoints[x_iter], y))  for x_iter in range(1, len(x_endpoints)) for y in y_endpoints)
        lines_list = tuple(l for l in itertools.chain(vertical_edges_iterator, horisontal_edges_iterator))
        for l in lines_list:
            l.set_growth_function(tessellation, symmetric_frame)
        for p in points_list:
            p.set_growth_function(tessellation, symmetric_frame)
        regions_lists.append(lines_list)
        regions_lists.append(points_list)

    return RegionsDescription(regions_lists, tessellation.dim, tessellation)