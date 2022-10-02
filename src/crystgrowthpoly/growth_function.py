from sage.all import *
import itertools
import functools
"""
This class represents growth function
"""
class growth_function:
    __PREDICATE_FORM_1_VARIABLE = "If n mod {} = {}, then"
    __PREDICATE_FORM_2_VARIABLES = "If k mod {} = {} and l mod {} = {}, then"
    def convert_coefficients_to_polynomial(self, coefficients,  number_of_variables):
        if number_of_variables != 1:
            indexes  = range(number_of_variables)
            row_scheme = tuple(itertools.chain(*(tuple(chosen_indexes for
                                                       chosen_indexes in itertools.combinations(indexes, number_of_variables - n))
                                                 for n in range(number_of_variables + 1))))
            return sum( coefficients[i] * functools.reduce(lambda acc, index,:  acc * self.__variables[index], row_scheme[i], 1)
                       for i in range(len(coefficients)))
        else:
            return sum(coefficients[-(1 + i)] * (n**i)  for i in range(len(coefficients)))

    def __init__(self, polynomials_coefficients, number_of_variables, cell_dimension, limits=None):
        if number_of_variables == 1:
            self.__variables = (var("n"))
        else:
            self.__variables = tuple(var('n_{}'.format(i)) for i in range (1, number_of_variables + 1))
        self.polynomials = tuple(
            self.convert_coefficients_to_polynomial(coef, number_of_variables)
            for coef in polynomials_coefficients)
        self.number_of_variables = number_of_variables
        self.cell_dimension = cell_dimension
        self.evaluation_function = lambda expr, arg: expr.subs(dict(zip(variables, arg)))
        if len(self.polynomials) == 1:
            self.is_polynomial = True
            self.polynomial_selector = lambda _: 0
        else:
            self.is_polynomial = False
            self.limits = tuple(limits)
            if number_of_variables == 1:
                self.N = limits
                self.polynomial_selector = lambda n: mod(n, limits[0])
            else:
                # TODO: In next commit rewrite this code to support higher dimmensions
                self.K = limits[0]
                self.L = limits[1]
                self.polynomial_selector = lambda arg: mod(arg[0], self.K) + N * mod(arg[1], self.L)
    def __call__(self, arg):
        return self.polynomials[self.polynomial_selector(arg)].subs(dict(zip(variables, arg)))
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
            # TODO: In next commit rewrite this code to support higher dimmensions
            var('k l')
            variables = tuple(var('n_{}'.format(i)) for i in range (1, self.number_of_variables + 1))
            if self.is_polynomial:
                pretty_print(fun, variables, "=", self.polynomials[0])
            else:
                for k_iter in range(self.K):
                    for l_iter in range(self.L):
                        print(growth_function.__PREDICATE_FORM_2_VARIABLES.format(self.K, k_iter, self.L, l_iter))
                        pretty_print(fun,(k,l), "=", self.polynomials[self.polynomial_selector((k_iter, l_iter))])