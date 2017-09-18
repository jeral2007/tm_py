import scipy as sc
#libcint params
PTR_LIGHT_SPEED    = 0
PTR_COMMON_ORIG    = 1
PTR_SHIELDING_ORIG = 4
PTR_RINV_ORIG      = 4
PTR_RINV_ZETA      = 7
PTR_ENV_START      = 20

CHARGE_OF  = 0
PTR_COORD  = 1
NUC_MOD_OF = 2
PTR_ZETA   = 3
RAD_GRIDS  = 4
ANG_GRIDS  = 5
ATM_SLOTS  = 6

ATOM_OF   = 0
ANG_OF    = 1
NPRIM_OF  = 2
NCTR_OF   = 3
KAPPA_OF  = 4
PTR_EXP   = 5
PTR_COEFF = 6
BAS_SLOTS = 8
def read_coord(filename):
    import scipy
    xyz, atoms = [], []
    for line in open(filename):
        if ('$coord' in line or line[0] == '#'):
            continue
        if '$end' in line:
            return xyz, atoms  # RETURN HERE
        aux = line.split()
        xyz += [scipy.array(map(float, aux[:3]))]
        atoms += [aux[-1]]
    raise ValueError('Incorrect coord file')


def read_embedding(filename, atoms):
    ii = 0
    qi, qq = [], []
    for line in open(filename):
        a, q = line.split()
        q = float(q)
        a = a.lower()
        while atoms[ii] != 'zz':
            ii +=1
        atoms[ii] = a
        qi += [ii]
        qq += [q]
    return qi, qq, atoms

class Geometry:
    def __init__(self, filename, filename_embedding=None):
        self.xyz, self.atoms = read_coord(filename)
        if filename_embedding is not None:
            self.qi, self.qq, self.atoms = read_embedding(filename_embedding, self.atoms)
    def atm(self, env, ptr=20):
        xyz = self.xyz
        res = sc.zeros((len(xyz), ATM_SLOTS), dtype=sc.int32)
        for i, r in enumerate(xyz):
            res[i,0:2] = [0, ptr+i*3]
            env[ptr + i*3:ptr+i*3+3] = r
        return res, ptr+len(xyz)*3

