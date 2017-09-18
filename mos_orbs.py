import read_basis as rb
import read_coord as rcd
import scipy as sc
"""explanation of symbols for d-,f-,and g-functions :
 d0  = (-xx-yy+2zz)/sqrt(12)
 d1a = xz
 d1b = yz
 d2a = xy
 d2b = (xx-yy)/2
 f0  = (-3xxx-3yyy+2zzz)/sqrt(60)
 f1a = (-xxx-xyy+4xzz)/sqrt(40)
 f1b = (-yyy-xxy+4yzz)/sqrt(40)
 f2a = xyz
 f2b = (xxz-yyz)/2
 f3a = (xxx-3xyy)/sqrt(24)
 f3b = (yyy-3xxy)/sqrt(24)
"""
fsymbols = {'s': ['s'], 'p': ['x', 'y', 'z'],
       'd': ['d0', 'd1a', 'd1b', 'd2a', 'd2b'],
       'f': ['f0', 'f1a', 'f1b', 'f2a', 'f2b',
             'f3a', 'f3b']}
def ljval(l):
    return "spdfgh".index(l)


def ao_order(basis, coord):

    ao_r = []
    ao_i = []
    ao_sym = []
    ao_alphas = []
    ao_coefs = []
    ao_spin  = []
    for spin in [-1, 1]:
        for a, r in zip(coord.atoms, coord.xyz):
            for i1 in basis.inds[a]:
                for fsym in fsymbols[basis.ls[i1]]:
                    ao_r += [r]
                    ao_i += [i1]
                    ao_sym += [fsym]
                    ao_alphas += [basis.alphas[i1]]
                    ao_coefs += [basis.coefs[i1]]
                    ao_spin  += [spin]
    return ao_r, ao_i, ao_sym, ao_alphas, ao_coefs, ao_spin

basis = rb.Basis('basis')
coord = rcd.Geometry('coord')
ii = 0
res = ao_order(basis, coord)
for r, i, sym, al, coe, sp in zip(*res):
    ii +=1
    print "{}. {} {} {}".format(ii, r, sym, sp)
env = sc.zeros(10000)
atm, ptr = coord.atm(env)
print atm
print ptr
print env[:ptr]
print "----------------------------"
bas, ptr = rb.make_bas(basis, coord, env, ptr)
sc.savetxt('bas.txt', fmt='%d', X=bas)
sc.savetxt('env.txt', env)
