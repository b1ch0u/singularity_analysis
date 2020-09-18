import itertools as it


def group_by_module(roots):
    '''
    Returns an iterator over groups of roots, grouped by increasing module.
    Assumes `roots` to be an iterable over 2-tuples of the form (root, multiplicity).
    '''
    def root_module(t):
        return abs(t[0])
    roots.sort(key=root_module)
    return it.groupby(roots, key=root_module)