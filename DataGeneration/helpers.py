import sympy as sym
import sympy.utilities.codegen as cgen
import sympy.codegen.rewriting as copt
import numpy as np
import re

# CREDITS TO HECTOR BARREIRO. Turns sympy code to C++ Eigen code.

def gradient(f, v):
    """
    Computes the gradient of f w.r.t. the given variables.
    """
    return sym.Matrix([f]).jacobian(v)

def hessian(f, v):
    """
    Computes the hessian of f w.r.t. the given variables.
    """
    return sym.hessian(f, v)

DTypeToEigen = { "double": "d", "float": "f", "int": "i" }
NDimsToEigen = [ "%s%.0s", "%.0sVectorX%s", "%.0sMatrixX%s" ]

def convert_var_to_eigen(var, lines, output=False):
    """
    Helper function to convert a generated SymPy variable to its corresponding
    Eigen type. Also performs the changes necessary to ensure column-major access
    to the data.
    """
    name = str(var.name)
    shape = [] if var.dimensions is None else [ 1+siz for _, siz in var.dimensions if siz != 0 ]
    ndims = len(shape)
    dtype = var.get_datatype('c')
    mtype = NDimsToEigen[ndims] % (dtype, DTypeToEigen[dtype])

    lines = [ line.replace(dtype + ' *' + str(var.name), ("const " if not output else "") + mtype + '& ' + str(var.name)) for line in lines ]

    if output is True and ndims > 0:
        index = np.where([ ("%s[0] =" % name) in line for line in lines ])[0][0]
        lines.insert(index, "    %s.resize(%s);" % (name, ", ".join([ str(size) for size in shape ])))

    if ndims == 2:
        for j in range(len(lines)):
            match = re.search('%s\[([0-9]+)\]' % name, lines[j])
            if match:
                idx = int(match.group(1))
                col = idx % shape[0]
                row = idx // shape[1]
                lines[j] = lines[j][:match.start(1)-1] + ("(%i, %i)" % (row, col)) + lines[j][match.end(1)+1:]

    return lines

def convert_if_else_to_inline(code):
    """
    Converts if-else and if-elseif-else structures to inline format.
    """
    pattern_if_elseif_else = r'const double if \((.*?)\) {\s*x(.*?) = (.*?);\s*}\s*elseif \((.*?)\) {\s*x\2 = (.*?);\s*}\s*else {\s*x\2 = (.*?);\s*}'
    pattern_if_else = r'const double if \((.*?)\) {\s*x(.*?) = (.*?);\s*}\s*else {\s*x\2 = (.*?);\s*}'

    replacement_if_elseif_else = r'const double x\2 = \1 ? \3 : (\4 ? \5 : \6);'
    replacement_if_else = r'const double x\2 = \1 ? \3 : \4;'

    code = re.sub(pattern_if_elseif_else, replacement_if_elseif_else, code, flags=re.DOTALL)
    code = re.sub(pattern_if_else, replacement_if_else, code, flags=re.DOTALL)

    return code

def generate(export, namestr):
    """
    Given a set of SymPy symbols, generates the corresponding C++ code to evaluate
    them in runtime.
    """
    generator = cgen.CCodeGen(cse=True)
    routines = []

    for i in range(len(export)):
        print(namestr+" ---    generate:: optimizing %i/%i ---" % (i+1, len(export)))
        name, expr = export[i]
        routine = generator.routine(name, expr, None, None)
        routines.append(routine)

    codes = []

    for i in range(len(routines)):
        print(" ---    generate:: writing %i/%i ---" % (i+1, len(routines)))
        routine = routines[i]
        code = generator.write((routine,), "", header=False, empty=True)
        codes.append(code)

    codes = [ code for (_, code), (_, _) in codes ]

    for i in range(len(codes)):
        lines = codes[i].split('\n')
        lines = lines[3:]
        lines = [ line.replace('   ', '    ') for line in lines ]
        for var in routines[i].arguments:
            lines = convert_var_to_eigen(var, lines, var in routines[i].result_variables)
        result = routines[i].result_variables[-1]
        name = str(result.name) if result.dimensions is not None else str(routines[i].name) + '_result'
        lines = [ line.replace(name, 'out') for line in lines ]
        codes[i] = '\n'.join(lines)

    return codes

def save(fname, str):
    """
    Saves the given string to a file.
    """
    with open(fname, "w") as f:
        f.write(str)
