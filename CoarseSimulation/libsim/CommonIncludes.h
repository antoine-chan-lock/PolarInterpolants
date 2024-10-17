#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <memory>
#include <vector>
#include <numeric>   // std::iota
#include <algorithm> // std::sort, std::stable_sort
#include <omp.h>

using namespace std;
using namespace Eigen;
using namespace Eigen::indexing;
