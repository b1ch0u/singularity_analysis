from sage.all import *


def get_lambda(k, i):
    return 1  # TODO


def get_e_k(k, alpha):
    return sum((-1)^i * get_lambda(k, i) * prod((alpha + 1 + i)
                                                for i in range(k))
               for i in range(k, 2*k + 1))


def gamma_derivatives(order, point, first_derivatives=[]):
    if not first_derivatives:
        s = var('s')
        first_derivatives.append(1 / gamma(-s))
    while order >= len(first_derivatives):
        first_derivatives.append(derivative(first_derivatives[-1]))
    return first_derivatives[order](point)