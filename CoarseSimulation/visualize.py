import numpy as np
import polyscope as ps

T = np.loadtxt("./mesh/T.csv", delimiter=',').astype(int)
X = np.loadtxt("./mesh/X.csv", delimiter=',')
Xdef = np.loadtxt("./mesh/Xdef.csv", delimiter=',')
ps.init()
ps.register_surface_mesh("rest", X, T)
ps.register_surface_mesh("def", Xdef, T)
ps.show()