from time import time

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

    ORDER = 4
    PRECISION = 1e-10

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
        except Exception as err:
            return 'Error: ' + str(err)

    def regular_tests():
        op = (1 + z) * Dz^2 + Dz
        print('\nlog(1+z) ->', my_extract(op, [0, 1, -1/2]))

        op = (1 - z) * Dz^2 - Dz
        print('\nlog(1-z) ->', my_extract(op, [0, -1, -1/2]))

        op = (z - 1) * Dz + 1
        print('\n1/(1-z) ->', my_extract(op, [1, 1]))

        op = (1 + z) * Dz + 1
        print('\n1/(1+z) ->', my_extract(op, [1, -1, 1, -1]))

        op = (1 - z^2) * Dz - 2*z
        print('\n1/(1 - z^2) ->', my_extract(op, [1, 0, 1, 0, 1]))

        op = (1 - z)^2 * Dz^2 - 3 * (1 - z) * Dz + 1
        print('\n1/(1-z) * log(1/(1-z)) (H_n) ->', my_extract(op, [0, 1, 3/2]))

        op = z * (1 - 2*z) * Dz - 1
        print('\nz/(1-2z) ->', my_extract(op, [0, 1, 2]))

        op = (1 + z^2) * Dz^2 + 2 * z * Dz
        print('\nArctan ->', my_extract(op, [0, 1, 0]))

        op = (4*z^2 - z)*Dz^2 + (14*z - 2)*Dz + 6
        print('\nrandom walks in Z*N ->', my_extract(op, [1, 3]))

        op = (16*z^4 - z^2)*Dz^3 + (128*z^3 + 8*z^2 - 6*z)*Dz^2 + (224*z^2 + 28*z - 6)*Dz + 64*z + 12
        print('\nrandom walks in N^2 ->', my_extract(op, [1,2,6,18]))

        op = (1 - z - z^2) * Dz^2 - (2 + 4*z) * Dz - 2
        print('\nFibonacci numbers ->', my_extract(op, [0, 1, 1, 2, 3, 5]))

        op = (4*z^2 - z)*Dz^2+(10*z-2)*Dz+2
        print('\nCatalan numbers ->', my_extract(op, [1, 1, 2, 5, 14]))

        op = (3*z^4 + 2*z^3 - z^2)*Dz^2 + (6*z^3 + 3*z^2 - z)*Dz + 1
        print('\nMotzkin numbers ->', my_extract(op, [1, 1, 2, 4, 9, 21, 51, 127], order=4))

    def apparent_sing_tests():
        print('\nTesting functions with an apparent singularity')
        op1 = 2 * (z - 1) * Dz - 3
        op2 = op1.lclm((z - 1/3) * Dz - 2)
        op3 = op1.lclm((z - 1/2) * Dz - 1)
        op4 = op2.lclm((z - 1/2) * Dz - 1)
        res1 = my_extract(op1, [1, -3/2, 3/8, 1/16])
        print('\n(1-z)^(3/2) ->', res1)
        res2 = my_extract(op2, [1, -3/2, 3/8, 1/16])
        print('\n(1-z)^(3/2) ->', res2)
        res3 = my_extract(op3, [1, -3/2, 3/8, 1/16])
        print('\n(1-z)^(3/2) ->', res3)
        res4 = my_extract(op4, [1, -3/2, 3/8, 1/16])
        print('\n(1-z)^(3/2) ->', res4)

        print('\nUndecidable singularity')
        op = (z^2 - 3*z + 2)*Dz^2 + (4*z - 6)*Dz + 2
        first_coefficients = [-5000000001/10^10, -2500000001/10^10, -1250000001/10^10]
        print('\n10^(-10)/(z-1) + 1/(z-2) (precision 1e-3) ->', my_extract(op, first_coefficients, precision=1e-5))
        print('\n10^(-10)/(z-1) + 1/(z-2) (precision 1e-20) ->', my_extract(op, first_coefficients, precision=1e-20))

    def no_sing_tests():
        print('\nTesting functions with no singularity, except maybe in 0.')

        op = (z^2-z)*Dz-(2*z-1)
        print('z(1-z) ->', my_extract(op, [0, 1, -1]))

        op = (-z^3+3*z^2-3*z+1)*Dz+(3*z^2-6*z+3)
        print('(z - 1)^3 ->', my_extract(op, [-1, 3, -3, 1]))

        op = (z - 1) * Dz - 1
        print('1-z ->', my_extract(op, [1, 1]))

        op = Dz - 1
        print('exp(z) ->', my_extract(op, [1, 1]))
        
        op = Dz^2 + 1
        print('sin(z) ->', my_extract(op, [0, 1, 0]))

    def irreg_tests():
        print('\nTesting functions with irregular singular points')

        op = z^2 * Dz + 1
        print('exp(1/z) ->', my_extract(op, [1,1]))

        op = (z - 1)^2 * Dz + 1
        print('exp(1/(z-1)) ->', my_extract(op, [1,1]))
    
    def weird_tests():
        print('\nCorner cases')
        print('(1-2z^3)Dz - 2z ->', my_extract((1-2*z^3)*Dz - 2*z, [1], verbose=False))

    def perf_tests():

        order = int(input('Order? '))
        precision = int(input('Precision? (10^___) (enter a negative number) '))
        N = int(input('Number of tests? '))

        op = (3*z^4 + 2*z^3 - z^2)*Dz^2 + (6*z^3 + 3*z^2 - z)*Dz + 1

        t1 = time()
        for _ in range(N): my_extract(op, [1, 1, 2, 4, 9, 21, 51, 127], order=20, precision=10^(precision))
        t2 = time()
        print(f'average time for Motzkin numbers up to order {order} and precision 10^{precision}: {(t2 - t1)/N}')

    regular_tests()
    apparent_sing_tests()
    no_sing_tests()
    irreg_tests()
    weird_tests()
    if input('\nRun performance test ? [y/N] ') == 'y':
        perf_tests()


if __name__ == '__main__':
    test_is_regular_singular_point()
    test_extract_asymptotics()