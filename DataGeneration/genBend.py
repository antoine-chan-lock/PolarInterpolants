import sympy as sym
import helpers as helpers

# This file contains a code generator for function, gradient and Hessian of the continuous beding operator.
# The symbolic output is generated in C++ in the file bend.inc.
# WARNING: sometimes in the generated code, .resize() is called within if(). It has tu be moved outside of the if() statement manually.

def generate(vparam):
    """
    Generate code for gradient, Hessian, and export routines.
    
    Args:
        vparam (list): List of parameters.
        
    Returns:
        str: Generated code.
    """
    U, vx, namestr = vparam[0], vparam[1], vparam[2]
    print(namestr + "--- Gradient generation ---")
    dUdF = helpers.gradient(U, [*vx])
    print(namestr + "--- Hessian generation ---")
    d2UdF2 = helpers.hessian(U, [*vx])
    
    # Define specific routines to export to code.
    export = [(namestr, U), ("d" + namestr + "dx", dUdF), ("d2" + namestr + "dx2", d2UdF2)]
    code = helpers.generate(export, namestr)
    return code


if __name__ == "__main__":
    x = sym.MatrixSymbol('x', 1, 3)
    s, c, a, k, d = sym.symbols('s c a k d')

    # In Plane
    S = sym.Matrix([[1 + s, 0, 0],
                    [0, 1 - c, 0],
                    [0, 0, 1]])
    U = sym.Matrix([[sym.cos(a), -sym.sin(a), 0],
                    [sym.sin(a), sym.cos(a), 0],
                    [0, 0, 1]])
    F = U * S * U.T

    xstretched = x * F

    # Bend
    R = sym.Matrix([[sym.cos(d), -sym.sin(d), 0],
                    [sym.sin(d), sym.cos(d), 0],
                    [0, 0, 1]])
    xbend = xstretched * R
    print(xbend.shape)
    xx = xbend[0, 0]
    yy = xbend[0, 1]
    zz = xbend[0, 2]

    kUpdated = k

    xcase1 = -sym.sin(kUpdated * xx) * (zz - 1 / kUpdated)
    zcase1 = sym.cos(kUpdated * xx) * (zz - 1 / kUpdated) + (1 / kUpdated)

    xcase2 = -xx * (kUpdated * zz - 1)
    zcase2 = zz + (kUpdated * (xx * xx)) / 2 - ((kUpdated * kUpdated) * (xx * xx) * zz) / 2

    # xnew = xcase1
    # znew = zcase1

    threshold = sym.symbols('threshold')
    cond = sym.Abs(kUpdated) > threshold
    xnew = sym.Piecewise((xcase1, cond), (xcase2, True))
    znew = sym.Piecewise((zcase1, cond), (zcase2, True))

    symx = generate([xnew, x, "symx"])
    symy = generate([yy, x, "symy"])
    symz = generate([znew, x, "symz"])

    concat = symx[0] + "\n" + symx[1] + "\n" + symx[2] + "\n" + symy[0] + "\n" + symy[1] + "\n" + symy[2] + "\n" + symz[0] + "\n" + symz[1] + "\n" + symz[2]

    concat = helpers.convert_if_else_to_inline(concat)

    helpers.save("./bend.inc", concat)
