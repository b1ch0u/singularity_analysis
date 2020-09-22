from ore_algebra import OreAlgebra

from ore_algebra.analytic.differential_operator import DifferentialOperator
from ore_algebra.analytic.local_solutions import (log_series, LocalBasisMapper,
                                                  simplify_exponent, LogMonomial)
from ore_algebra.analytic.path import Point

from utils import group_by_module


DEFAULT_PRECISION = 1e-53
DEFAULT_ORDER = 5


def is_regular_singular_point(point, leading_mult, z, op):
    '''
    Checks weather a point is a regular singular point of a differential operator.
    `leading_mult` is the multiplicity of `point` as a root
    of the leading coefficient of `op`.
    '''
    r = op.order()
    for i, poly in enumerate(op):
        if r - i < leading_mult \
                and not ((z - point)^(leading_mult - (r - i))).divides(poly):
            return False
    return True


def my_expansions(op, point, order):
    mypoint = Point(point, op)
    ldop = op.shift(mypoint)
    class Mapper(LocalBasisMapper):
        def fun(self, ini):
            return log_series(ini, self.shifted_bwrec, order)
    sols = Mapper(ldop).run()
    return [[(c / ZZ(k).factorial(),
              point,
              sol.leftmost + n,
              k)
             for n, vec in enumerate(sol.value)
             for k, c in reversed(list(enumerate(vec)))
             if c != 0]
            for sol in sols]


def compute_initial_decomp(op, first_coefficients):
    '''
    Compute the decomposition of the solution to differential operator `op`
    with first coefficients `first_coefficients` in the local basis in 0,
    as calculated by local_basis_expansions.
    '''
    distinguished_monomials = Point(0, op).local_basis_structure()
    proj = [0] * op.order()
    for i, sol in enumerate(distinguished_monomials):
        n, k = simplify_exponent(sol.valuation), sol.log_power
        if n >= 0 and k == 0:
            try:
                proj[i] = first_coefficients[n]
            except IndexError:
                raise Exception('Not enough first coefficients supplied')
    return proj


def handle_monomial(zeta, alpha, k, order, verbose, result_var):
    '''
    Compute an asymptotic expansion of the coefficients of
    (z-zeta)^alpha * (log (z-zeta))^k
    up to `order` terms.
    '''
    corr = 2 * pi * I  #  * e^(alpha * corr)
    if k == 0:
        res = (-zeta)^alpha * simplify(asymptotic_expansions.SingularityAnalysis(result_var,
                                                zeta=zeta,
                                                alpha=-alpha,
                                                precision=order,
                                                normalized=False))
        if verbose:
            print(f'1. got called with {zeta}, {alpha}, {k} and returned ', res)
        return res
    # print((log(-1/zeta))^(1-1))
    # tmp_res = asymptotic_expansions.SingularityAnalysis(result_var,
    #                                             zeta=zeta,
    #                                             alpha=-alpha,
    #                                             beta=k,
    #                                             precision=order,
    #                                             normalized=False)
    # print('tmp res', tmp_res)
    # print(dir(tmp_res))
    # print(simplify_exponent(tmp_res))
    # print(simplify(tmp_res))
    res = (-1)^k * (
            sum(binomial(k, i) * (log(-1/zeta))^(k-i) * simplify(asymptotic_expansions.SingularityAnalysis(result_var,
                zeta=zeta,
                alpha=-alpha,
                beta=i,
                precision=order,
                normalized=False)) for i in range(1, k+1)))
    if alpha not in NN:
        res += (log(-1/zeta))^k
    res *= (-zeta)^alpha
    if verbose:
        print(f'2. got called with {zeta}, {alpha}, {k} to order {order} and returned ', res)
    return res


def handle_root(op, root, ini, order, precision, verbose, result_var):
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

    local_expansions = my_expansions(op, root, order=order + 3)  # TODO eviter de refaire le calcul deja fait dans transition_matrix
    if verbose:
        print('basis in', root, ':', op.local_basis_expansions(root, order))
        print('coeffs in local basis', coeffs_in_local_basis)

    res = []
    found_a_sing = False
    for lambda_i, expansion in zip(coeffs_in_local_basis, local_expansions):
        for i, (c, _, alpha, m) in enumerate(expansion[:order]):
            if alpha not in NN or m > 0:
                found_a_sing = True
                res.append(lambda_i * c * handle_monomial(root,
                                                       alpha,
                                                       m,
                                                       order - i + 1,
                                                       verbose,
                                                       result_var))
            else:
                res.append(0)
        if order < len(expansion):
            _, _, alpha, m = expansion[order]
            last_term_expansion = handle_monomial(root, alpha, m, 0, verbose, result_var)
            res.append(last_term_expansion)
    if verbose:
        print('contribution in', root, ':', sum(res))
    return res, found_a_sing


def extract_asymptotics(op,
                        first_coefficients,
                        order=DEFAULT_ORDER,
                        precision=DEFAULT_PRECISION,
                        verbose=False,
                        result_var='n'):
    '''
    Compute an asymptotic expansion of a solution to the differential operator `op`
    specified by its first coefficients `first_coefficients`, with
    constants certified up to arbitrary precision, and up to any desired order.

    Assuming `z` is the variable used in `op`.
    `order` controls the order of the expansion.
    `precision` controls the disired amount of certified precision for the constants.
    '''
    op = DifferentialOperator(op)
    z = op.base_ring().gen()
    ini = compute_initial_decomp(op, first_coefficients)
    if verbose:
        print('decomp in 0:', ini)
        print('basis in 0:', op.local_basis_expansions(0))
    leading_roots = op.leading_coefficient().roots(QQbar)
    # TODO check 0 is reg or sing reg ?

    all_monomials = []
    can_stop = False
    for module, roots in group_by_module(leading_roots):
        if module == 0:
            continue

        for root, multiplicity in roots:
            if not is_regular_singular_point(root, multiplicity, z, op):
                    raise ValueError('Got an irregular singular point. Stopping here')
            contribution, found_a_sing = handle_root(op, root, ini, order, precision, verbose, result_var)
            can_stop = can_stop or found_a_sing
            all_monomials.extend(contribution)

        if can_stop:
            break
    if all_monomials:
        return sum(all_monomials)
    raise ValueError('No singularity was found')