#include "Mesh.h"
#include "Tools.h"
#include "CommonIncludes.h"

// If on macos, run: export KMP_DUPLICATE_LIB_OK=TRU
// then : ./build/example 1 3
// then : python visualize.py
// 1 stands for material 1, 3 stands for in-plane stretching boundary condition

void runSingle(int idmat, int idBC)
{
    string dir = "./mesh/";
    Mesh mesh(dir, idmat, idBC);
    mesh.initElements();
    mesh.setMaterialInit(1);
    mesh.clampDoFs();
    MatrixXd X = mesh.getRest();
    MatrixXd x(X.rows(), X.cols() + 1);
    x << X, MatrixXd::Zero(X.rows(), 1);
    VectorXd vX = Tools::vectorize(x);
    bool success;
    vX = mesh.newtonRaphson(vX, success);
    mesh.setMaterialInit(0);
    vX = mesh.newtonRaphson(vX, success);
    Tools::writeMatrix(Tools::devectorize(vX), dir + "/Xdef.csv");
}

int main(int argc, char *argv[])
{
    int idmat = atof(argv[1]);
    int idBC = atof(argv[2]);
    cout << "idmat : " << idmat << endl;
    cout << "idBC : " << idBC << endl;
    runSingle(idmat, idBC);
    return 0;
}
