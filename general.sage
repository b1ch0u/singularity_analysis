from ore_algebra import OreAlgebra
from ore_algebra.analytic.differential_operator import DifferentialOperator
from ore_algebra.analytic.local_solutions import (log_series, LocalBasisMapper,
                                                  simplify_exponent, LogMonomial)
from ore_algebra.analytic.path import Point

from utils import group_by_module


DEFAULT_PRECISION = 1e-53
DEFAULT_ORDER = 5

SA = asymptotic_expansions.SingularityAnalysis


def assert_point_is_regular_singular(point, leading_mult, op):
    '''
    Checks weather a point is a regular singular point of a differential operator.
    `leading_mult` is the multiplicity of `point` as a root
    of the leading coefficient of `op`.
    '''
    z = op.base_ring().gen()
    r = op.order()
    for i, poly in enumerate(op):
        if r - i < leading_mult \
                and not ((z - point)^(leading_mult - (r - i))).divides(poly):
            raise ValueError(f'Got an irregular singular point ({point}).')


def check_0_is_regular_singular_point(op, leading_roots):
    if any(root == 0 for root, _ in leading_roots):
        mult_of_zero = next(mult for root, mult in leading_roots if root == 0)
        assert_point_is_regular_singular(0, mult_of_zero, op)


def my_expansions(op, point, order):
    '''
    Compute the local basis expansion of `op` in the neighbourhood of `point`.
    Same job as the "real" local_basis_expansions functions, but returns
    a list of coefficients and powers, instead of a FormalSum, since these
    are not easily manipulated.
    '''
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
        if n in NN and n >= 0 and k == 0:
            try:
                proj[i] = first_coefficients[n]
            except IndexError:
                raise Exception('Not enough first coefficients supplied')
    return proj


def handle_monomial(zeta, alpha, k, order, verbose, result_var, precision):
    '''
    Compute an asymptotic expansion of the coefficients of
    (z-zeta)^alpha * (log (z-zeta))^k
    up to `order` terms.
    '''
    if k == 0:
        res = (-zeta)^alpha * SA(result_var,
                                 zeta=zeta,
                                 alpha=-alpha,
                                 precision=order,
                                 normalized=False)
        if verbose:
            print(f'1. got called with zeta={zeta}, alpha={alpha}, and returned ', res)
        return res

    res = (-1)^k * sum(binomial(k, i) * (log(-1/zeta))^(k-i) * SA(result_var,
                                                                  zeta=zeta,
                                                                  alpha=-alpha,
                                                                  beta=i,
                                                                  precision=order,
                                                                  normalized=False)
                        for i in range(1, k+1))
    if alpha not in NN:
        res += (log(-1/zeta))^k * handle_monomial(zeta, alpha, 0, order, verbose, result_var, precision)
    res *= (-zeta)^alpha
    if verbose:
        print(f'2. got called with zeta={zeta}, alpha={alpha}, k={k} to order {order} and returned ', res)
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

    local_expansions = my_expansions(op, root, order=order + 3)
    if verbose:
        print('basis in', root, ':', op.local_basis_expansions(root, order + 3))
        print('coeffs in local basis', coeffs_in_local_basis)
        print('local expansions', local_expansions)

    res = []
    found_a_sing = False
    for lambda_i, expansion in zip(coeffs_in_local_basis, local_expansions):
        for i, (c, _, alpha, m) in enumerate(expansion[:order]):
            if alpha not in NN or m > 0:

                if alpha not in QQ:
                    raise ValueError('Got an irrational value for alpha')

                if lambda_i.is_nonzero():
                    found_a_sing = True
                res.append(lambda_i * c * handle_monomial(root,
                                                          alpha,
                                                          m,
                                                          order - i + 1,
                                                          verbose,
                                                          result_var,
                                                          precision))
            else:
                res.append(0)

        if order < len(expansion):
            _, _, alpha, m = expansion[order]
            if alpha not in NN or m > 0:
                last_term_expansion = handle_monomial(root, alpha, m, 0, verbose, result_var, precision)
                res.append(last_term_expansion)
    
    if verbose: print('contribution in', root, ':', sum(res))
    
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

    leading_roots = op.leading_coefficient().roots(QQbar)

    check_0_is_regular_singular_point(op, leading_roots)

    ini = compute_initial_decomp(op, first_coefficients)
    if verbose:
        print('basis in 0:', op.local_basis_expansions(0))
        print('decomp in 0:', ini)

    all_monomials_by_modulus = []
    can_stop = False
    for module, roots in group_by_module(leading_roots):
        if module == 0:
            continue

        monomials = []
        for root, multiplicity in roots:
            assert_point_is_regular_singular(root, multiplicity, op)

            contribution, found_a_sing = handle_root(op, root, ini, order, precision, verbose, result_var)
            can_stop = can_stop or found_a_sing
            monomials.extend(contribution)

        if verbose: print('found a sing:', found_a_sing)

        all_monomials_by_modulus.append(sum(monomials))

        if can_stop:
            break

    if any(monomials for monomials in all_monomials_by_modulus):
        return (all_monomials_by_modulus)
    raise ValueError('No singularity was found')