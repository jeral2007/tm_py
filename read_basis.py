import scipy as sc
from math import *
def _read_alphas(stream, ngtos):
    a, c = sc.zeros(ngtos), sc.zeros(ngtos)
    for n in xrange(ngtos):
        line = stream.next()
        [a[n], c[n]] = map(float, line.split())
    if ngtos ==1:
        c[0] = 1
    return a,c


def read_basis(filename):
    state = 0 #wait for $basis
    atom_is, ls, alphas, coefs = [], [], [], []
    fl = open(filename)
    for line in fl:
        if "$ecp" in line:
            break #end of basis definition
        if line[0] == '#':
            continue # skip comments
        if state == 0:
            if "$basis" in line:
                state = 1 # wait for asterisk
                continue
        elif state == 1:
            if "*" in line:
                state = 2 # read atom name
                continue
            else:
                raise ValueError("Incorrect basis file {}, line {}".format(filename, line))
        elif state == 2:
            atom_name, _ = line.split()
            state = 3 # wait for asterisk
            continue
        elif state == 3:
            if "*" in line:
                state = 4 # read basis function
                continue
            else:
                raise ValueError("Incorrect basis file {}, line {}".format(filename, line))
        elif state == 4:
            if '*' in line:
                state=2 #next atom
                continue
            ngto, l = line.split()
            atom_is += [atom_name]
            ls += [l]
            aa, cc = _read_alphas(fl, int(ngto)) # using fact, that each value in iterator is returned only once
            alphas += [aa]
            coefs += [cc]
    return atom_is, ls, alphas, coefs

class Basis(object):
    def __init__(self, filename):
        self.atom_is, self.ls, self.alphas, self.coefs = read_basis(filename)
        self.inds = {}
        for i, a in enumerate(self.atom_is):
            if a in self.inds:
                self.inds[a] += [i]
            else:
                self.inds[a] = [i]

def gto_norm(n, a):
    # normalization factor of function r^n e^{-a r^2}
    s = 2**(2*n+3) * factorial(n+1) * (2*a)**(n+1.5) \
            / (factorial(2*n+2) * sqrt(pi))
    return sc.sqrt(s)


def make_bas(basis, coord, env, ptr):
    bas = []
    cptr = ptr
    subarr = lambda a, inds: [a[i] for i in inds]
    for i, (r, a) in enumerate(zip(coord.xyz, coord.atoms)):
        inds = basis.inds[a]
        ls = map("spdfgh".index, subarr(basis.ls, inds))
        alphas = subarr(basis.alphas, inds)
        coefs = subarr(basis.coefs, inds)
        for l in xrange(min(ls), max(ls)+1):
            for ii in xrange(len(ls)):
                if ls[ii] != l: continue
                al = alphas[ii]
                cc = coefs[ii]
                bas +=[[i, l, len(al), 1, 0, cptr, cptr+len(al), 0]]
                env[cptr:cptr+len(al)] = al
                cptr += len(al)
                env[cptr:cptr+len(cc)] = gto_norm(l, al) * cc
                cptr += len(al)
    return sc.array(bas, dtype=sc.int32), cptr
