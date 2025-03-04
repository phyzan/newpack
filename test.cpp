#include "ode.hpp"

using Tt = double;
using Tf = vec<Tt, 4>;

Tf f(const Tt& t, const Tf& q, const std::vector<Tt>& args){
    return {q[2], q[3], -q[0], -q[1]};
    // return q;
}

Tt fevent(const Tt& t, const Tf& f, const std::vector<Tt>& args){
    return f[1]-1;
}

bool check_fevent(const Tt& t, const Tf& f, const std::vector<Tt>& args){
    return f[3]>0;
}


int main(){

    double pi = 3.14159265359;
    double t_max = 10001*pi/2;
    Tf q0 = {1, 1, 2.3, 4.5};

    // SolverArgs<Tt, Tf, true, true> S = {f, 0., 1000, q0, 1e-3, 0., 1e-8, 0., {}, fevent, nullptr, check_fevent, nullptr, nullptr, nullptr, nullptr, 1e-12};

    ODE<Tt, 4> ode(f, 0, q0, 1e-2, 1e-5, 1e-10, 0., {}, "RK23", 1e-10, nullptr, nullptr, nullptr);
    // ode.integrate(t_max/2, -1, 20, false);
    // ode.integrate(t_max/2, -1, 20, false);
    // OdeResultReference<Tt, Tf> res = ode.integrate(t_max, 10, 5, false);
    ode.integrate(t_max, 2);
    ode.state().show();
    std::cout << ode.runtime << "\n";
    // for (size_t i=1; i<ode.t.size(); i++){
    //     std::cout << ode.t[i]-ode.t[i-1] << "\n";
    // }
    // while (true){
    //     ode.advance().show();
    //     std::cin.get();
    // }

    // res.examine();

    // while (true){
        
    // }

    // res.examine();
    // ode.state().show();

    // for (size_t i=0; i<ode.t.size(); i++){
    //     std::cout << ode.t[i] << "\n";
    // }

    // for (size_t i=0; i<res.events.size(); i++){
    //     std::cout << res.q[res.events[i]].transpose() << "\n";
    // }
    // std::cout << ode.result().q[200023];

}
// #include <functional>


// template <class Tt, class Ty>
// class Function{

// private:

//     Ty(*_func_ptr)(const Tt&, const Ty&) = nullptr;
//     std::function<Ty(const Tt&, const Ty&)> _std_func;


// public:

//     Function(std::function<Ty(const Tt&, const Ty&)> f){
//         _std_func = f;
//     }

//     Function(Ty(*f)(const Tt&, const Ty&)){
//         _func_ptr = f;
//     }

//     Ty operator()(const Tt& t, const Ty& y){
//         return (_func_ptr != nullptr) ? _func_ptr(t, y) : _std_func(t, y);
//     }
// };





// #include <eigen3/Eigen/Dense>
// #include <iostream>

// int main() {
//     const int rows = 5;
//     const int cols = 4;

//     // Create an Eigen matrix and array
//     Eigen::Matrix<double, rows, cols> mat;
//     Eigen::Array<double, 1, cols> arr;

//     // Fill the array with values (you can use setRandom() for random values)
//     arr.setRandom();

//     // Print the original array
//     std::cout << "Array:\n" << arr << std::endl;

//     // Assign the first row of the matrix to the array
//     mat.row(0) = arr.matrix();  // Convert the array to a matrix and assign it

//     // Print the result
//     std::cout << "Matrix after assignment:\n" << mat << std::endl;

//     return 0;
// }
