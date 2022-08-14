from sage.all import *
"""
This class represents growth function
"""
class growth_function:
    __PREDICATE_FORM_1_VARIABLE = "If n mod {} = {}, then"
    __PREDICATE_FORM_2_VARIABLES = "If k mod {} = {} and l mod {} = {}, then"
    def convert_coefficients_to_polynomial(coefficients,  number_of_variables):
        if number_of_variables == 2:
            var('k l')
            return coefficients[3] + l * coefficients[2] + k * coefficients[1] + (k*l) * coefficients[0]
        else:
            var('n')
            return coefficients[2] + n * coefficients[1] + (n**2) * coefficients[0]

    def __init__(self, polynomials_coefficients, number_of_variables, cell_dimension, limits=None):
        self.polynomials = tuple(
            growth_function.convert_coefficients_to_polynomial(coef, number_of_variables) for coef in polynomials_coefficients)
        self.number_of_variables = number_of_variables
        self.cell_dimension = cell_dimension
        if number_of_variables == 1:
            var("n")
            self.evaluation_function = lambda expr, n_val: expr.subs({n: n_val})
        else:
            var("k l")
            self.evaluation_function = lambda expr, arg: expr.subs({k: arg[0], l: arg[1]})
        if len(self.polynomials) == 1:
            self.is_polynomial = True
        else:
            self.is_polynomial = False
            if number_of_variables == 1:
                self.N = limits
                self.polynomial_selector = lambda n: mod(n, limits)
            else:
                self.K = limits[0]
                self.L = limits[1]
                self.polynomial_selector = lambda arg: mod(arg[0], self.K) + N * mod(arg[1], self.L)
    def __call__(self, arg):
        if self.is_polynomial:
            return self.evaluation_function(self.polynomials[0], (arg))
        else:
            return self.evaluation_function(self.polynomials[self.polynomial_selector(arg)], arg)
    def show(self):
        fun = var("f_{}".format(self.cell_dimension))
        if self.number_of_variables == 1:
            var('n')
            if self.is_polynomial:
                pretty_print(fun,"(",n,")", "=", self.polynomials[0])
            else:
                for n_iter in range(self.N):
                    print(growth_function.__PREDICATE_FORM_1_VARIABLE.format(self.N, n_iter))
                    pretty_print(fun,"(",n,")", "=", self.polynomials[n_iter])
        else:
            var('k l')
            if self.is_polynomial:
                pretty_print(fun,(k,l), "=", self.polynomials[0])
            else:
                for k_iter in range(self.K):
                    for l_iter in range(self.L):
                        print(growth_function.__PREDICATE_FORM_2_VARIABLES.format(self.K, k_iter, self.L, l_iter))
                        pretty_print(fun,(k,l), "=", self.polynomials[self.polynomial_selector((k_iter, l_iter))])