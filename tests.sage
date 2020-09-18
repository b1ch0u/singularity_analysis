from sage.all import *

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


if __name__ == '__main__':
    test_is_regular_singular_point()

    #####

    Pols.<z> = PolynomialRing(QQ)
    Diff.<Dz> = OreAlgebra(Pols)

    order = 4
    precision = 1e-10

    def my_extract(op, first_coefficients):
        return extract_asymptotics(op, first_coefficients, z, order=order, precision=precision)

    op = (1 - z - z^2) * Dz^2 - (2 + 4*z) * Dz - 2
    print('Fibonacci numbers ->', my_extract(op, [0, 1, 1, 2, 3, 5]))

    op = (z - 1) * Dz^2 + Dz
    print('log(1-z) ->', my_extract(op, [0, -1, -1/2]))

    op = (1 - z^2) * Dz - 2*z
    print('1/(1 - z^2) ->', my_extract(op, [1, 0, 1, 0, 1]))

    op = (1 + z^2) * Dz^2 + 2 * z * Dz
    print('Arctan', my_extract(op, [0, 1, 0])) 

    op = (4*z^2 - z)*Dz^2 + (14*z - 2)*Dz + 6
    print('random walks in a half plane :', my_extract(op, [1, 3]))

    op = (16*z^4 - z^2)*Dz^3 + (128*z^3 + 8*z^2 - 6*z)*Dz^2 + (224*z^2 + 28*z - 6)*Dz + 64*z + 12
    print('random walks in quadrant :', my_extract(op, [1,2,6,18]))

    op = (1 - z) * Dz - 1
    print('1/(1-z) ->', my_extract(op, [1, 1]))

    op = (1 + z) * Dz + 1
    print('1/(1+z) ->', my_extract(op, [1, -1, 1, -1]))

    op = -2/3 * (1 - z) * Dz - 1
    print('(1-z)^(3/2) ->', my_extract(op, [1, 1]))

    op = z * (1 - 2*z) * Dz - 1
    print('z/(1-2z) ->', my_extract(op, [0, 1, 2]))

    print('\nTesting functions with no singularity, except maybe in 0.\n')

    op = Dz - 1
    try:
        print('exp(z) ->', my_extract(op, [1, 1]))
    except ValueError as e:
        print('exp(z) ->', str(e))

    op = (z - 1) * Dz - 1
    try:
        print('1-z ->', my_extract(op, [1, 1]))
    except ValueError as e:
        print('1-z ->', str(e))
    
    op = Dz^2 + 1
    try:
        print('sin(z) ->', my_extract(op, [0, 1, 0]))
    except ValueError as e:
        print('sin(z) ->', str(e))
    
    op = z * Dz^2 - 2 * z^2 * Dz - 4*z
    try:
        print('Catalan numbers ->', my_extract(op, [1, 1, 2, 5, 14]))
    except ValueError as e:
        print('Catalan numbers ->', str(e))