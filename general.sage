#from sage.all import *

from ore_algebra import *

from ore_algebra.analytic.differential_operator import DifferentialOperator
from ore_algebra.analytic.local_solutions import (log_series, LocalBasisMapper,
                                                  simplify_exponent, LogMonomial)
from ore_algebra.analytic.path import Point

import itertools as it

from utils import *


EPS = 1e-20
ORDER = 4


def group_by_module(roots):
    '''
    Returns an iterator over groups of roots, grouped by increasing module.
    Assumes `roots` to be an iterable over 2-tuples of the form (root, multiplicity).
    '''
    def root_module(t):
        return abs(t[0])
    roots.sort(key=root_module)
    return it.groupby(roots, key=root_module)


def is_regular_singular_point(point, leading_mult, var, op):
    '''
    Checks wether a point is a regular singular point of a differential operator.
    `leading_mult` is the multiplicity of `point` as a root
    of the leading coefficient of `op`.
    '''
    r = op.order()
    for i, poly in enumerate(op): # 0 -> p_0, ...
        if r - i < leading_mult \
                and not ((var - point)^(leading_mult - (r - i))).divides(poly):
            return False
    return True


def my_expansions(op, point, order=None, ring=None):
    mypoint = Point(point, op)
    dop = DifferentialOperator(op)
    ldop = dop.shift(mypoint)
    if order is None:
        ind = ldop.indicial_polynomial(ldop.base_ring().gen())
        order = max(dop.order(), ind.dispersion()) + 3
    class Mapper(LocalBasisMapper):
        def fun(self, ini):
            return log_series(ini, self.shifted_bwrec, order)
    sols = Mapper(ldop).run()
    return [[(c / ZZ(k).factorial() * (CBF(-point))^(CBF(sol.leftmost + n)),
                point,
                sol.leftmost + n,
                k)
                for n, vec in enumerate(sol.value)
                for k, c in reversed(list(enumerate(vec)))]
            for sol in sols]


class GeneralMonomial:
    def __init__(self, c, zeta, alpha, m):
        # TODO verifier domaines d'analyticite
        self.c = c
        self.zeta = zeta
        self.alpha = alpha
        self.m = m
        self.asymp_c = c / gamma(-alpha)

    def is_singular(self):
        return self.alpha < 0 or self.m != 0

    def __repr__(self):
        if self.m == 0:
            return f'{self.asymp_c} * {1/self.zeta}^n * n^{-self.alpha - 1}'
        return f'{self.asymp_c} * {1/self.zeta}^n * n^{-self.alpha - 1} * log^{self.m}(n)'

    def asymptotic_symbolic(self):
        n = var('n')
        return (self.c / gamma(-self.alpha)) * self.zeta^n * n^(-self.alpha - 1) * log(n)^self.m

    def __cmp__(self, other):
        biggest = None
        if abs(self.zeta) != abs(other.zeta):
            biggest = self if abs(self.zeta) < abs(other.zeta) else other  # le plus petit |zeta| gagne
        elif self.alpha != other.alpha:
            biggest = self if self.alpha < other.alpha else other  # le plus petit alpha gagne
        elif self.m != other.m:
            biggest = self if self.m > other.m else other  # le plus grand m gagne
        if biggest is None or not biggest.c.is_nonzero():
            return 0
        return 1 if biggest is self else -1

    def __le__(self, other):
        return 0 if self.__cmp__(other) == 1 else 1

    def __lt__(self, other):
        return 1 if self.__cmp__(other) == -1 else 0


class Monomial:
    def __init__(self, coeff, exp_base, n_power, log_power):
        self.coeff = coeff
        self.exp_base = exp_base
        self.n_power = n_power
        self.log_power = log_power


def basic_transfer(zeta, alpha, order):
    res = []
    common_coeff = (-zeta)^alpha / gamma(-alpha)
    for k in range(order):
        res.append(Monomial(common_coeff * get_e_k(k, alpha),
                            1/zeta,
                            -alpha - 1 - k,
                            0))
    return res


def log_transfer(zeta, alpha, m, order):
    # TODO renormalisation
    pass


def compute_term_expansion(coeff, alpha, m, zeta, order):
    if alpha in NN and m == 0:
        return []
    if alpha not in NN and m == 0:
        return basic_transfer(zeta, alpha, order)
    return log_transfer(zeta, alpha, m, order)


# class Expansion:
#     def __init__(self):
#         self.terms = 0
#         self.bound = 0

#     def add_term(self, new_term, term_bound):
#         self.terms += new_term
#         self.bound += term_bound

#     def __repr__(self):
#         return f'{self.terms} +-{self.bound}'


def compute_initial_decomp(op, first_coefficients):
    def my_local_monomials(op, point):
        dop = DifferentialOperator(op)
        struct = Point(point, dop).local_basis_structure()
        x = SR(dop.base_ring().gen()) - point
        return [(
                    1 / sol.log_power.factorial(),
                    simplify_exponent(sol.valuation),
                    sol.log_power)
                for sol in struct]

    proj = [0] * op.order()
    for i, (alpha, n, k) in enumerate(my_local_monomials(op, 0)):
        if k == 0 and n >= 0:
            proj[i] = first_coefficients[n]
    return proj


def handle_root(op, root, ini):
    trans_matrix = op.numerical_transition_matrix([0, root],
                                                  eps=EPS,
                                                  assume_analytic=True)
    coeffs_in_local_basis = trans_matrix * vector(ini)

    local_expansions = my_expansions(op, root, order=ORDER)

    new_monomials = [GeneralMonomial(lambda_i * c, zeta, alpha, m)
                     for terms in local_expansions
                     for lambda_i, (c, zeta, alpha, m) in zip(coeffs_in_local_basis, terms)
                     if alpha not in NN or m]  # exclude the polynomial case
    if not new_monomials:
        print('got no new monomial in', root)

    return new_monomials


def extract_asymptotics(op, first_coefficients):
    ini = compute_initial_decomp(op, first_coefficients)

    leading_roots = op.leading_coefficient().roots(QQbar)
    # TODO check 0 is reg or sing reg

    found_a_sing = False
    all_monomials = []
    for module, roots in group_by_module(leading_roots):
        roots = list(roots)
        if module == 0:
            continue

        for root, multiplicity in roots:

            if not is_regular_singular_point(root, multiplicity, z, op):
                    raise ValueError('Got an irregular singular point. Stopping here')

            new_monomials = handle_root(op, root, ini)
            all_monomials += new_monomials

            if any(mon.is_singular() and mon.c.is_nonzero()
                   for mon in new_monomials):
                found_a_sing = True

        if found_a_sing:
            all_monomials.sort(reverse=True)
            return sum((mon.asymptotic_symbolic() for mon in all_monomials), start=0)


def test_is_regular_singular_point():
    Pols.<z> = PolynomialRing(QQ)
    Diff.<Dz> = OreAlgebra(Pols)

    op = (4*z^2 - z) * Dz^2 + (14 * z - 2) * Dz + 6
    assert is_regular_singular_point(0, 1, z, op)
    assert is_regular_singular_point(1/4, 1, z, op)

    op2 = (4*z^3 - z^2) * Dz^2 + (14 * z - 2) * Dz + 6
    assert not is_regular_singular_point(0, 2, z, op2)
    assert is_regular_singular_point(1/4, 1, z, op2)

    op3 = (z)^4 * (z-1)^5 * Dz^3 + (z^2) * (z - 1) * Dz
    assert is_regular_singular_point(0, 4, z, op3)
    assert not is_regular_singular_point(1, 5, z, op3)

    op4 = (z-1/4)^8 * (z-1/5)^12 * Dz^9 + (z - 1/4)^4 * Dz^5 + 1
    assert is_regular_singular_point(1/4, 8, z, op4)
    assert not is_regular_singular_point(1/5, 12, z, op4)


if __name__ == '__main__':
    test_is_regular_singular_point()

    #####

    Pols.<z> = PolynomialRing(QQ)
    Diff.<Dz> = OreAlgebra(Pols)

    # random walks in a half plane
    op = (4*z^2 - z)*Dz^2 + (14*z - 2)*Dz + 6
    print(extract_asymptotics(op, [1, 3]))

    # random walks in quadrant
    op = (16*z^4 - z^2)*Dz^3 + (128*z^3 + 8*z^2 - 6*z)*Dz^2 + (224*z^2 + 28*z - 6)*Dz + 64*z + 12
    print(extract_asymptotics(op, [1,2,6,18]))