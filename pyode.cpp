#include "pyode.hpp"

PYBIND11_MODULE(_integrate, m) {
    define_ode_module<double, Eigen::Array<double, 1, Eigen::Dynamic>>(m);
}
