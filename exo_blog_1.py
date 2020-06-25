from functools import lru_cache


@lru_cache(maxsize=None)
def count_walks_xy(x, y, n):
    '''
    Count walks that stay in the upper-half plane.
    '''
    if y < 0:
        return 0
    if n == 0 and x == 0 and y == 0:
        return 1
    if n == 0:
        return 0
    return count_walks_xy(x - 1, y, n - 1) + \
           count_walks_xy(x + 1, y, n - 1) + \
           count_walks_xy(x, y - 1, n - 1) + \
           count_walks_xy(x, y + 1, n - 1)


def count_walks(n):
    return sum([count_walks_xy(x, y, n)
                for x in range(-n, n + 1)
                for y in range(n + 1)])


LST = [count_walks(k) for k in range(11)]
print(LST)

'''
OE = oeis(LST)[0]
print(OE)
print(OE.comments())
'''

import time

t1 = time.time()
w_100 = count_walks(100)
t_tot = time.time() - t1
print('w_100 :', w_100)
print('computed in', t_tot, 'seconds')


t1 = time.time()
w_150 = count_walks(150)
t_tot = time.time() - t1
print('w_150 :', w_150)
print('computed in', t_tot, 'seconds')


'''
from ore_algebra import *

# Create the shift ore_algebra with rational coefficients
Ind.<n> = PolynomialRing(QQ)
Shift.<Sn> = OreAlgebra(Ind)

# Guess a recurrence satisfied by the sequence (need about 30 terms to find the recurrence)
LST = [count_walks(k) for k in range(30)]
rec = guess(LST, Shift)
print(rec)

LongLST = rec.to_list(LST,301)
print('w_300 :', LongLST[300])

print('prob of staying in the upper half plane for 300 steps :', (LongLST[300]/4^(300)).n())


# Compute the first two asymptotic terms of a basis of the recurrence
print('Une base de l\'espace des solutions de la rec')
print(rec.generalized_series_solutions(2))


print('Une etude heuristique pour approximer lambda_1')
LongLST = rec.to_list(LST,10001)
print((LongLST[100]/(4^100/100)).n())
print((LongLST[1000]/(4^1000/1000)).n())
print((LongLST[10000]/(4^10000/10000)).n())



### La partie qui me concerne commence vraiment maintenant


# Create the differential algebra to encode linear differential equation
Pols.<z> = PolynomialRing(QQ)
Diff.<Dz> = OreAlgebra(Pols)

# Guess an annihilating linear differential equation for generating function
diffWalk = guess(LST,Diff)
print(diffWalk)


# Converting from the differential equation to a recurrence 
# gives rec (up to a left multiple)
print(Diff(diffWalk).to_S(Shift))
'''
