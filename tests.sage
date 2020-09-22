from sage.all import *  # TODO

import os
os.system('sage --preparse general.sage')
os.system('mv general.sage.py general.py')
from general import *

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


def test_extract_asymptotics():
    Pols.<z> = PolynomialRing(QQ)
    Diff.<Dz> = OreAlgebra(Pols)

    ORDER = 5
    PRECISION = 1e-5

    def my_extract(op, first_coefficients,
                       verbose=False, order=ORDER, precision=PRECISION,
                       result_var='n'):
        try:
            return extract_asymptotics(op,
                                   first_coefficients,
                                   order=order,
                                   precision=precision,
                                   verbose=verbose,
                                   result_var=result_var)
        except ValueError as e:
            return str(e)

    # op = (z - 1) * Dz^2 + Dz
    # #print('\nlog(1-z) ->', my_extract(op, [0, -1, -1/2], verbose=False))
    # res = my_extract(op, [0, -1, -1/2])
    # print(res.symbolic_expression())
    # print(res.exact_part())
    # #print(dir(res.exact_part().symbolic_expression()))
    # n = var('n')
    # wanted_exact = -n^(-1) + n^(-2) - n^(-3) + n^(-4) - n^(-5)
    # print(dir(wanted_exact))
    # exact = res.exact_part().symbolic_expression()
    # print(wanted_exact)
    # print('coeff n -2', wanted_exact.coefficient(n^(-2)))
    # print('approx n -2', res.symbolic_expression().coefficient(n^(-2)))
    # print('1 in cbf', 1 in CBF(res.symbolic_expression().coefficient(n^(-2))))
    # print(wanted_exact.operands()[0].coefficient(1/n))
    # print(wanted_exact.coefficients(n))

    # op = (1 + z) * Dz^2 + Dz
    # print('\nlog(1+z) ->', my_extract(op, [0, 1, -1/2], verbose=False, order=4))

    # op = (1 - z) * Dz^2 - Dz
    # print('\nlog(1-z) ->', my_extract(op, [0, -1, -1/2], verbose=False, order=4))

    # op = (z - 1) * Dz + 1
    # print('\n1/(1-z) ->', my_extract(op, [1, 1]))

    # op = (1 + z) * Dz + 1
    # print('\n1/(1+z) ->', my_extract(op, [1, -1, 1, -1]))

    # op = (1 - z^2) * Dz - 2*z
    # print('\n1/(1 - z^2) ->', my_extract(op, [1, 0, 1, 0, 1]))

    # op = (1 - z)^2 * Dz^2 - 3 * (1 - z) * Dz + 1
    # print('\n1/(1-z) * log(1/(1-z)) (H_n) ->', my_extract(op, [0, 1, 3/2], verbose=False, order=2, result_var='N'))

    # op = z * (1 - 2*z) * Dz - 1
    # print('\nz/(1-2z) ->', my_extract(op, [0, 1, 2]))

    # op = (1 + z^2) * Dz^2 + 2 * z * Dz
    # print('\nArctan', my_extract(op, [0, 1, 0]))

    # op = (4*z^2 - z)*Dz^2 + (14*z - 2)*Dz + 6
    # print('\nrandom walks in Z*N :', my_extract(op, [1, 3]))

    # op = (16*z^4 - z^2)*Dz^3 + (128*z^3 + 8*z^2 - 6*z)*Dz^2 + (224*z^2 + 28*z - 6)*Dz + 64*z + 12
    # print('\nrandom walks in N^2 :', my_extract(op, [1,2,6,18], order=7, verbose=False))

    # op = (1 - z - z^2) * Dz^2 - (2 + 4*z) * Dz - 2
    # print('\nFibonacci numbers ->', my_extract(op, [0, 1, 1, 2, 3, 5], verbose=False))

    # op = (4*z^2 - z)*Dz^2+(10*z-2)*Dz+2
    # print('\nCatalan numbers ->', my_extract(op, [1, 1, 2, 5, 14], verbose=False))

    # op = (3*z^4 + 2*z^3 - z^2)*Dz^2 + (6*z^3 + 3*z^2 - z)*Dz + 1
    # print('\nMotzkin numbers ->', my_extract(op, [1, 1, 2, 4, 9, 21, 51, 127], verbose=False, order=4))

    # only introducing an apparent singularity
    op1 = 2 * (z - 1) * Dz - 3
    op2 = op1.lclm(z * Dz - 3)
    op3 = op1.lclm((z - 1/2) * Dz - 1)
    op4 = op2.lclm((z - 1/2) * Dz - 1)
    res1 = my_extract(op1, [1, -3/2, 3/8, 1/16], verbose=True)
    res2 = my_extract(op2, [1, -3/2, 3/8, 1/16], verbose=True)
    res3 = my_extract(op3, [1, -3/2, 3/8, 1/16], verbose=True)
    res4 = my_extract(op4, [1, -3/2, 3/8, 1/16], verbose=True)
    print('\n(1-z)^(3/2) ->', res1)
    print('\n(1-z)^(3/2) ->', res2)
    print('\n(1-z)^(3/2) ->', res3)
    print('\n(1-z)^(3/2) ->', res4)

    #####
    print('\nTesting functions with no singularity, except maybe in 0.\n')

    op = (z^2-z)*Dz-(2*z-1)
    print('\nz(1-z) ->', my_extract(op, [0, 1, -1], verbose=False))

    op = (-z^3+3*z^2-3*z+1)*Dz+(3*z^2-6*z+3)
    print('\n(z - 1)^3 ->', my_extract(op, [-1, 3, -3, 1]))

    op = (z - 1) * Dz - 1
    print('\n1-z ->', my_extract(op, [1, 1]))

    op = Dz - 1
    print('\nexp(z) ->', my_extract(op, [1, 1]))
    
    op = Dz^2 + 1
    print('\nsin(z) ->', my_extract(op, [0, 1, 0]))


if __name__ == '__main__':
    #test_is_regular_singular_point()
    test_extract_asymptotics()