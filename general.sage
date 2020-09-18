from ore_algebra import *

from ore_algebra.analytic.differential_operator import DifferentialOperator
from ore_algebra.analytic.local_solutions import (log_series, LocalBasisMapper,
                                                  simplify_exponent, LogMonomial)
from ore_algebra.analytic.path import Point

from utils import group_by_module


DEFAULT_PRECISION = 1e-10
DEFAULT_ORDER = 5


def is_regular_singular_point(point, leading_mult, z, op):
    '''
    Checks wether a point is a regular singular point of a differential operator.
    `leading_mult` is the multiplicity of `point` as a root
    of the leading coefficient of `op`.
    '''
    r = op.order()
    for i, poly in enumerate(op):
        if r - i < leading_mult \
                and not ((z - point)^(leading_mult - (r - i))).divides(poly):
            return False
    return True


def my_expansions(op, point, order=None):
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
    return [[(c / ZZ(k).factorial(),
                point,
                sol.leftmost + n,
                k)
                for n, vec in enumerate(sol.value)
                for k, c in reversed(list(enumerate(vec)))]  # TODO if c != 0
            for sol in sols]


def compute_initial_decomp(op, first_coefficients):
    '''
    Compute the decomposition of the solution to differential operator `op`
    with first coefficients `first_coefficients` in the local basis in 0,
    as calculated by local_basis_expansions.
    '''
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


def handle_monomial(zeta, alpha, k, order):
    '''
    Compute an asymptotic expansion of the coefficients of
    (z-zeta)^alpha * (log (z-zeta))^k
    up to `order` terms.
    '''
    if k == 0:
        return (-zeta)^alpha * asymptotic_expansions.SingularityAnalysis('n',
                                                zeta=zeta,
                                                alpha=-alpha,
                                                precision=order)
    return (-zeta)^alpha * sum(binomial(k, i) * log(-1/zeta)^(k-i) * asymptotic_expansions.SingularityAnalysis('n',
                                                zeta=zeta,
                                                alpha=-alpha,
                                                beta=i,
                                                precision=order) for i in range(1, k+1))


def handle_root(op, root, ini, order, precision):
    '''
    Compute the contribution to asymptotic expansion of the coefficients f_n
    of a solution f to a given differential operator `op`, caracterized by
    a list `ini` of first coefficients in expansion at 0,
    in the neighbourhood of `root`.    
    '''
    trans_matrix = op.numerical_transition_matrix([0, root],
                                                  eps=precision,
                                                  assume_analytic=True)
    coeffs_in_local_basis = trans_matrix * vector(ini)

    local_expansions = my_expansions(op, root, order=order)  # TODO eviter de refaire le calcul deja fait dans transition_matrix

    return [lambda_i * c * handle_monomial(root, alpha, m, order)  # reduce order when not necessary
                     for lambda_i, expansion in zip(coeffs_in_local_basis, local_expansions)
                     for (c, _, alpha, m) in expansion
                     if c != 0 and (alpha not in NN or m > 0)]


def extract_asymptotics(op,
                        first_coefficients,
                        z,
                        order=DEFAULT_ORDER,
                        precision=DEFAULT_PRECISION):
    '''
    Compute an asymptotic expansion of a solution to the differential operator `op`
    specified by its first coefficients `first_coefficients`, with
    constants certified up to arbitrary precision, and up to any desired order.

    Assuming `z` is the variable used in `op`.
    `order` controls the order of the expansion.
    `precision` controls the disired amount of certified precision for the constants.
    '''
    ini = compute_initial_decomp(op, first_coefficients)
    leading_roots = op.leading_coefficient().roots(QQbar)
    # TODO check 0 is reg or sing reg ?

    all_monomials = []
    for module, roots in group_by_module(leading_roots):
        if module == 0:
            continue

        for root, multiplicity in roots:
            if not is_regular_singular_point(root, multiplicity, z, op):
                    raise ValueError('Got an irregular singular point. Stopping here')
            all_monomials += handle_root(op, root, ini, order, precision)

        if all_monomials:
            return sum(all_monomials)
    raise ValueError('No singularity was found')