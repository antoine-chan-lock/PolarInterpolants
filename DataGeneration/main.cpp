
#include "DataGeneration.h"

int main(int argc, char *argv[])
{
   int nmat = atof(argv[1]);
   cout << nmat << endl;

   DataGeneration dg(to_string(nmat));
   dg.sampling();
   cout << "Sampling done" << endl;
   dg.genDerivatives();
   cout << "Derivatives done" << endl;
}
