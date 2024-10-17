import sympy as sym
import helpers as helpers

def generate(vparam):
    U,vx,namestr = vparam[0], vparam[1], vparam[2]
    dUdF,d2UdF2 = helpers.derivativesParallel(U, [ *vx ])
    export = [(namestr, U), ("d"+namestr+"dx", dUdF), ("d2"+namestr+"dx2", d2UdF2)]
    code = helpers.generate(export,namestr)
    helpers.save("symStrains/"+namestr+".inc", helpers.convert_if_else_to_inline(code[0]))
    helpers.save("symStrains/d"+namestr+"dx.inc", helpers.convert_if_else_to_inline(code[1]))
    helpers.save("symStrains/d2"+namestr+"dx2.inc", helpers.convert_if_else_to_inline(code[2]))

def generatef(vparam):
    U,vx,namestr = vparam[0], vparam[1], vparam[2]
    export = [(namestr, U)]
    code = helpers.generate(export,namestr)
    helpers.save("symStrains/"+namestr+".inc", helpers.convert_if_else_to_inline(code[0]))

def norm(symx):
    return sym.sqrt( symx.dot(symx) )



if __name__ == "__main__":
    vx = sym.MatrixSymbol('vx', 18, 1) #x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6 z1 z2 z3 z4 z5 z6
                                       #0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17
    vXrest = sym.MatrixSymbol('vXrest', 6, 1) #X1 X2 X3 Y1 Y2 Y3
    vexist = sym.MatrixSymbol('vexist', 3, 1)

    x = sym.Matrix([[vx[0], vx[3], vx[6]],
                    [vx[1], vx[4], vx[7]],
                    [vx[2], vx[5], vx[8]]])
    xflap = sym.Matrix([[vx[9], vx[12], vx[15]],
                        [vx[10], vx[13], vx[16]],
                        [vx[11], vx[14], vx[17]]])
    Xrest = sym.Matrix([[vXrest[0], vXrest[2], vXrest[4]],
                        [vXrest[1], vXrest[3], vXrest[5]]])


    #In Plane strains
    Dm = sym.Matrix.hstack(Xrest[:,1] - Xrest[:,0], Xrest[:,2] - Xrest[:,0])
    Ds = sym.Matrix.hstack(x[:,1] - x[:,0], x[:,2] - x[:,0])
    F = Ds * Dm.inv()
    E2 = (F.T * F - sym.eye(2))
    e1,e2,a = helpers.eig(E2)
    s = sym.sqrt( e1 + 1 ) -1
    c = -sym.sqrt( e2 + 1 ) +1

    # Bending strains
    K = sym.Matrix.zeros(2,2)
    e01 = x[:,1] - x[:,0]
    e02 = x[:,2] - x[:,0]
    n0 = e01.cross(e02)
    R2_col1 = F[:,0]
    R2_col1 = R2_col1 / norm(R2_col1)
    R2_col2 = n0.cross(F[:,0])
    R2_col2 = R2_col2 / norm(R2_col2)
    R2 = sym.Matrix.hstack(R2_col1, R2_col2)
    R2 = R2.T


    for i in range(3):
        xa = x[:,i]
        xb = x[:,(i+1)%3]
        xn = xflap[:,i]
        t1 = xn-xa
        t2 = xb-xa
        ni = t1.cross(t2)
        ei = xa - xb
        theta = helpers.dihedral(n0,ni,ei)
        ti = R2 * ( (-n0).cross(ei) )
        M = (ti * ti.T) / ( norm(ei)*norm(n0)*n0.dot(n0) )
        K += vexist[i]*theta*M

    k1,k2,d = helpers.eig(K)

    from multiprocessing import Process
    vparam = [[s,vx,"sym_s"],
              [c,vx,"sym_c"],
              [a,vx,"sym_a"],
              [k1,vx,"sym_k1"],
              [k2,vx,"sym_k2"],
              [d,vx,"sym_d"]]

    # generatef() to generate the function, generate() to generate the function + gradient + hessian

    generatef(vparam[0])
    generatef(vparam[1])
    generatef(vparam[2])
    generatef(vparam[3])
    generatef(vparam[4])
    generatef(vparam[5])