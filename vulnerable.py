# Written in Python 3
# exec(open("./edet.py").read()) 
import numpy as np
import scipy as sp
import scipy.sparse as sparse
import scipy.sparse.csgraph as csgraph
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
infile = '../Preprocess/Sets21corr/21sn15_c40'
outfile = './vulnerable.out'
with open(infile) as f:
  lines = (line for line in f if not line.startswith('#'))
  U0, V0 = np.loadtxt(lines, delimiter='\t', unpack=True, dtype=int)

U = U0
V = V0
n = max(np.maximum(U, V)) + 1
e = U.size
A =  sparse.csr_matrix((np.ones(e), (U, V)), shape=(n, n), dtype=bool)
A += sparse.csr_matrix((np.ones(e), (V, U)), shape=(n, n), dtype=bool)
k = np.squeeze(np.asarray(A.sum(axis=0))) # np.matrix -> np.array
k_large = np.maximum(k[U], k[V])
k_max = max(k)
k_hist = np.histogram(k, bins=np.arange(k_max+2))[0]
nk = np.count_nonzero(k_hist)
print("No. nodes:", n, ", No. edges", e, ", Max. deg: ", k_max, ", Non-0 deg: ", nk)

n_vul = np.zeros(nk, dtype=int)
e_vul = np.zeros(nk, dtype=int)
k_hist_nz = np.nonzero(k_hist)[0]
for i in range(nk):
  phiinv = k_hist_nz[i]
  U = U0[k_large <= phiinv];
  V = V0[k_large <= phiinv];
  e = U.size
  A =  sparse.csr_matrix((np.ones(e), (U, V)), shape=(n, n), dtype=bool)
  A += sparse.csr_matrix((np.ones(e), (V, U)), shape=(n, n), dtype=bool)

  n_comp, vec_comp = csgraph.connected_components(A, directed=False)
  comp_hist = np.histogram(vec_comp, bins=np.arange(n_comp))[0]
  giant_ind = np.argmax(comp_hist)
  n_vul[i] = comp_hist[giant_ind]

  gidx = np.arange(n)[vec_comp==giant_ind]
  G = A.tocsr()[gidx, :].tocsc()[:, gidx]
  giant_list = sparse.find(G)
  e_vul[i] = giant_list[0].size/2

outfile = open(outfile,'w')
print("[", end='', file=outfile)
for i in range(nk):
  print("[", k_hist_nz[i], ", ", n_vul[i], ", ", e_vul[i], "]", sep='', end='', file=outfile)
  if i < nk - 1:
    print(", ", end='', file=outfile)

print("]", end='', file=outfile)
outfile.close()
