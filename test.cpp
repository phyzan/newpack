#include "ode.hpp"

using Tt = double;
using Tf = vec<Tt, 4>;

Tf f(const Tt& t, const Tf& q, const std::vector<Tt>& args){
    return {q[2], q[3], -q[0], -q[1]};
    // return q;
}

Tt fevent(const Tt& t, const Tf& f, const std::vector<Tt>& args){
    return f[1]-2;
}

Tt fevent2(const Tt& t, const Tf& f, const std::vector<Tt>& args){
    return f[1]-3;
}

bool check_fevent(const Tt& t, const Tf& f, const std::vector<Tt>& args){
    return f[3]>0;
}

Tf mask(const Tt& t, const Tf& f, const std::vector<Tt>& args){
    return {1, 1, 0, 0};
}


int main(){

    double pi = 3.14159265359;
    double t_max = 10001*pi/2;
    Tf q0 = {1, 1, 2.3, 4.5};

    // SolverArgs<Tt, Tf, true, true> S = {f, 0., 1000, q0, 1e-3, 0., 1e-8, 0., {}, fevent, nullptr, check_fevent, nullptr, nullptr, nullptr, nullptr, 1e-12};

    Event<Tt, Tf> event1("Event1", fevent, nullptr, mask);
    StopEvent<Tt, Tf> event2("Event2", fevent);

    ODE<Tt, Tf> ode(f, 0, q0, 1e-2, 1e-6, 1e-12, 1e-8, {}, "RK45", 1e-10, {event1});
    // ode.integrate(t_max/2, -1, 20, false);
    // ode.integrate(t_max/2, -1, 20, false);
    // OdeResultReference<Tt, Tf> res = ode.integrate(t_max, 10, 5, false);
    // OdeResult<Tt, Tf> res1 = ode.integrate(t_max/2, -1, 5).as_copy();
    // OdeResult<Tt, Tf> res2 = ode.integrate(t_max/2, -1, 5).as_copy();
    // res2.examine();

    // for (const auto& [ev, array] : res1.events){
    //     std::cout << std::endl;
    //     for (const size_t& i : array){
    //         std::cout << res1.t[i] << " q: " << res1.q[i] << "\n";
    //     }
    // }
    // for (const auto& [ev, array] : res2.events){
    //     std::cout << std::endl;
    //     for (const size_t& i : array){
    //         std::cout << res2.t[i] << " q: " << res2.q[i] << "\n";
    //     }
    // }
    // std::cout << ode.runtime << "\n";
    // ode.state().show();
    // for (size_t i=1; i<ode.t.size(); i++){
    //     std::cout << ode.t[i]-ode.t[i-1] << "\n";
    // }
    ode.integrate(100).examine();
    while (true){
        ode.advance().show();
        std::cin.get();
    }

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

// #ifndef PYODE_HPP
// #define PYODE_HPP


// #include <pybind11/pybind11.h>
// #include <pybind11/numpy.h>
// #include <functional>
// #include "ode.hpp"

// namespace py = pybind11;

// template<class Tt, class Ty>
// Func<Tt, Ty> to_Func(const py::object f, const std::vector<size_t>& shape);

// std::vector<size_t> shape(const py::array& arr) {
//     const ssize_t* shape_ptr = arr.shape();  // Pointer to shape data
//     size_t ndim = arr.ndim();  // Number of dimensions
//     return std::vector<size_t>(shape_ptr, shape_ptr + ndim);
// }

// template<class Tt, class Ty>
// event_f<Tt, Ty> to_event(const py::object py_event, const std::vector<size_t>& shape);

// template<class Tt, class Ty>
// is_event_f<Tt, Ty> to_event_check(const py::object py_event_check, const std::vector<size_t>& shape);

// py::dict to_PyDict(const std::map<std::string, std::vector<size_t>>& _map) {
//     py::dict py_dict;
//     for (const auto& [key, vec] : _map) {
//         py::array_t<size_t> np_array(vec.size(), vec.data()); // Create NumPy array
//         py_dict[key.c_str()] = np_array; // Assign to dictionary
//     }
//     return py_dict;
// }

// template<class Tt, class Ty>
// Ty toCPP_Array(const py::array& A);


// template<class Scalar, class ArrayType>
// py::array_t<Scalar> to_numpy(const ArrayType& array, const std::vector<size_t>& _shape = {});

// template<class Tt, class Ty>
// std::vector<Tt> flatten(const std::vector<Ty>&);

// template<class T>
// py::tuple to_tuple(const std::vector<T>& vec);

// template<class Tt, class Ty>
// std::vector<Event<Tt, Ty>> to_Events(const py::object& events, const std::vector<size_t>& shape);


// template<class Tt, class Ty>
// std::vector<StopEvent<Tt, Ty>> to_StopEvents(const py::object& events, const std::vector<size_t>& shape);


// #pragma GCC visibility push(hidden)

// template<class Tt, class Ty>
// class PyEvent : public Event<Tt, Ty>{

// public:

//     PyEvent(py::str name,
//             py::object when,
//             py::object check_if,
//             py::object mask,
//             const std::vector<size_t>& shape={})
//             :Event<Tt, Ty>(name.cast<std::string>(), to_event<Tt, Ty>(when, shape), to_event_check<Tt, Ty>(check_if, shape), to_Func<Tt, Ty>(mask, shape)), _py_when(when), _py_check_if(check_if), _py_mask(mask){}

//     PyEvent<Tt, Ty> remake(const std::vector<size_t>& shape){
//         return PyEvent<Tt, Ty>(this->name, _py_when, _py_check_if, _py_mask, shape);
//     }

// private:
    
//     py::object _py_when;
//     py::object _py_check_if;
//     py::object _py_mask;
// };

// template<class Tt, class Ty>
// class PyStopEvent : public StopEvent<Tt, Ty>{

// public:

//     PyStopEvent(py::str name,
//             py::object when,
//             py::object check_if,
//             const std::vector<size_t>& shape={})
//             :StopEvent<Tt, Ty>(name.cast<std::string>(), to_event<Tt, Ty>(when, shape), to_event_check<Tt, Ty>(check_if, shape)), _py_when(when), _py_check_if(check_if) {}

//     PyStopEvent<Tt, Ty> remake(const std::vector<size_t>& shape){
//         return PyStopEvent<Tt, Ty>(this->name, _py_when, _py_check_if, shape);
//     }

// private:
    
//     py::object _py_when;
//     py::object _py_check_if;
// };


// template<class Tt, class Ty>
// class PySolverState : public SolverState<Tt, Ty>{

// public:
//     PySolverState(const Tt& t, const Ty& q, const Tt& habs, const std::string& event, const bool& diverges, const bool& is_stiff, const bool& is_running, const bool& is_dead, const size_t& N, const std::string& message, const std::vector<size_t>& shape): SolverState<Tt, Ty>(t, q, habs, event, diverges, is_stiff, is_running, is_dead, N, message), shape(shape) {}

//     const std::vector<size_t> shape;

// };


// template<class Tt, class Ty>
// class PyODE : public ODE<Tt, Ty>{
// public:

//     PyODE(py::object f, const Tt t0, const py::array q0, const Tt stepsize, const Tt rtol, const Tt atol, const Tt min_step, const py::tuple args, const py::str method, const Tt event_tol, py::object events, py::object stop_events) : ODE<Tt, Ty>(
//         to_Func<Tt, Ty>(f, shape(q0)), t0, toCPP_Array<Tt, Ty>(q0), stepsize, rtol, atol, min_step, 
//         toCPP_Array<Tt, std::vector<Tt>>(args), 
//         method.cast<std::string>(), event_tol, to_Events<Tt, Ty>(events, shape(q0)), to_StopEvents<Tt, Ty>(stop_events, shape(q0))), _shape(shape(q0)) {}

//     PySolverState<Tt, Ty> py_state() const{
//         SolverState<Tt, Ty> s = this->state();
//         return PySolverState<Tt, Ty>(s.t, s.q, s.habs, s.event, s.diverges, s.is_stiff, s.is_running, s.is_dead, s.N, s.message, _shape); 
//     }

//     std::vector<size_t> full_shape() const{
//         std::vector<size_t> result;
//         std::vector<size_t> _shape = shape(this->q[0]);
//         result.reserve(1 + _shape.size()); // Pre-allocate memory for efficiency
//         result.push_back(this->t.size());        // Add the first element
//         result.insert(result.end(), _shape.begin(), _shape.end()); // Append the original vector
//         return result;
//     }

// private:
//     const std::vector<size_t> _shape;
// };



// #pragma GCC visibility pop


// /*
// -----------------------------------------------------------------------
// -----------------------------IMPLEMENTATIONS-------------------------------
// -----------------------------------------------------------------------
// */


// template<class Tt, class Ty>
// Func<Tt, Ty> to_Func(const py::object f, const std::vector<size_t>& shape){
//     if (f.is_none()){
//         return nullptr;
//     }
//     Func<Tt, Ty> g = [f, shape](const Tt& t, const Ty& y, const std::vector<Tt>& args) -> Ty {
//         return toCPP_Array<Tt, Ty>(f(t, to_numpy<Tt>(y, shape), *to_tuple(args)));
//     };
//     return g;
// }

// template<class Tt, class Ty>
// event_f<Tt, Ty> to_event(const py::object py_event, const std::vector<size_t>& shape){
//     if (py_event.is_none()){
//         return nullptr;
//     }
//     event_f<Tt, Ty> g = [py_event, shape](const Tt& t, const Ty& f, const std::vector<Tt>& args) -> Tt {
//         return py_event(t, to_numpy<Tt>(f, shape), *to_tuple(args)).template cast<Tt>();
//     };
//     return g;
// }

// template<class Tt, class Ty>
// is_event_f<Tt, Ty> to_event_check(const py::object py_event_check, const std::vector<size_t>& shape){
//     if (py_event_check.is_none()){
//         return nullptr;
//     }
//     is_event_f<Tt, Ty> g = [py_event_check, shape](const Tt& t, const Ty& f, const std::vector<Tt>& args) -> bool {
//         return py_event_check(t, to_numpy<Tt>(f, shape), *to_tuple(args)).equal(py::bool_(true));
//     };
//     return g;
// }

// template<class Tt, class Ty>
// Ty toCPP_Array(const py::array& A){
//     size_t n = A.size();
//     Ty res(n);

//     const Tt* data = static_cast<const Tt*>(A.data());

//     for (size_t i=0; i<n; i++){
//         res[i] = data[i];
//     }
//     return res;
// }

// template<class Scalar, class ArrayType>
// py::array_t<Scalar> to_numpy(const ArrayType& array, const std::vector<size_t>& _shape){
//     if (_shape.size() == 0){
//         return py::array_t<Scalar>(shape(array), array.data());
//     }
//     else{
//         return py::array_t<Scalar>(_shape, array.data());
//     }
// }

// template<class Tt, class Ty>
// std::vector<Tt> flatten(const std::vector<Ty>& f){
//     size_t nt = f.size();
//     size_t nd = f[0].size();
//     std::vector<Tt> res(nt*nd);

//     for (size_t i=0; i<nt; i++){
//         for (size_t j=0; j<nd; j++){
//             res[i*nd + j] = f[i][j];
//         }
//     }
//     return res;
// }

// template<class T>
// py::tuple to_tuple(const std::vector<T>& arr) {
//     py::tuple py_tuple(arr.size());  // Create a tuple of the same size as the vector
//     for (size_t i = 0; i < arr.size(); ++i) {
//         py_tuple[i] = py::float_(arr[i]);  // Convert each double element to py::float_
//     }
//     return py_tuple;
// }





// template<class Tt, class Ty>
// void define_ode_module(py::module& m) {
//     py::class_<PyEvent<Tt, Ty>>(m, "Event", py::module_local())
//         .def(py::init<py::str, py::object, py::object, py::object>(),
//             py::arg("name"),
//             py::arg("when"),
//             py::arg("check_if")=py::none(),
//             py::arg("mask")=py::none());

//     py::class_<PyStopEvent<Tt, Ty>>(m, "StopEvent", py::module_local())
//         .def(py::init<py::str, py::object, py::object>(),
//             py::arg("name"),
//             py::arg("when"),
//             py::arg("check_if")=py::none());

//     py::class_<OdeResult<Tt, Ty>>(m, "OdeResult", py::module_local())
//         .def_property_readonly("t", [](const OdeResult<Tt, Ty>& self){
//             return to_numpy<Tt>(self.t);
//         })
//         .def_property_readonly("q", [](const OdeResult<Tt, Ty>& self){
//             return to_numpy<Tt>(flatten<Tt, Ty>(self.q), self.full_shape());
//         })
//         .def_property_readonly("events", [](const OdeResult<Tt, Ty>& self){
//             return to_PyDict(self.events);
//         })
//         .def_property_readonly("diverges", [](const OdeResult<Tt, Ty>& self){return self.diverges;})
//         .def_property_readonly("is_stiff", [](const OdeResult<Tt, Ty>& self){return self.is_stiff;})
//         .def_property_readonly("success", [](const OdeResult<Tt, Ty>& self){return self.success;})
//         .def_property_readonly("runtime", [](const OdeResult<Tt, Ty>& self){return self.runtime;})
//         .def_property_readonly("message", [](const OdeResult<Tt, Ty>& self){return self.message;})
//         .def("examine", &OdeResult<Tt, Ty>::examine);


//     py::class_<PySolverState<Tt, Ty>>(m, "SolverState", py::module_local())
//         .def_property_readonly("t", [](const PySolverState<Tt, Ty>& self){return self.t;})
//         .def_property_readonly("q", [](const PySolverState<Tt, Ty>& self){return to_numpy<Tt>(self.q, self.shape);})
//         .def_property_readonly("event", [](const PySolverState<Tt, Ty>& self){return (self.event == "") ? py::none() : self.event;})
//         .def_property_readonly("diverges", [](const PySolverState<Tt, Ty>& self){return self.diverges;})
//         .def_property_readonly("is_stiff", [](const PySolverState<Tt, Ty>& self){return self.is_stiff;})
//         .def_property_readonly("is_running", [](const PySolverState<Tt, Ty>& self){return self.is_running;})
//         .def_property_readonly("is_dead", [](const PySolverState<Tt, Ty>& self){return self.is_dead;})
//         .def_property_readonly("N", [](const PySolverState<Tt, Ty>& self){return self.N;})
//         .def_property_readonly("message", [](const PySolverState<Tt, Ty>& self){return self.message;})
//         .def("show", &PySolverState<Tt, Ty>::show);

        

//     py::class_<PyODE<Tt, Ty>>(m, "LowLevelODE", py::module_local())
//         .def(py::init<py::object, Tt, py::array, Tt, Tt, Tt, Tt, py::tuple, py::str, Tt, py::object, py::object, py::object, py::object, py::object, py::object, py::object>(),
//             py::arg("f"),
//             py::arg("t0"),
//             py::arg("q0"),
//             py::arg("stepsize"),
//             py::kw_only(),
//             py::arg("rtol")=1e-6,
//             py::arg("atol")=1e-12,
//             py::arg("min_step")=0.,
//             py::arg("args")=py::tuple(),
//             py::arg("method")="RK45",
//             py::arg("event_tol")=1e-12,
//             py::arg("events")=py::none(),
//             py::arg("stop_event")=py::none())
//         .def("integrate", [](PyODE<Tt, Ty>& self, const Tt& interval, const int& max_frames, const int& max_events, const bool& terminate, const bool& display){
//             return self.integrate(interval, max_frames, max_events, terminate, display);
//             },
//             py::arg("interval"),
//             py::kw_only(),
//             py::arg("max_frames")=-1,
//             py::arg("max_events")=-1,
//             py::arg("terminate")=true,
//             py::arg("display")=false)
//         .def("advance", &PyODE<Tt, Ty>::advance)
//         .def("state", &PyODE<Tt, Ty>::py_state)
//         .def_property_readonly("t", [](const PyODE<Tt, Ty>& self){return to_numpy<Tt>(self.t);})
//         .def_property_readonly("q", [](const PyODE<Tt, Ty>& self){return flatten<Tt, Ty>(self.q), self.full_shape();})
//         .def("events", [](const PyODE<Tt, Ty>& self){return to_PyDict(self.event_map());})
//         .def_property_readonly("runtime", [](const PyODE<Tt, Ty>& self){return self.runtime;});

// }


// //g++ -O3 -Wall -shared -std=c++20 -fopenmp -I/usr/include/python3.12 -I/usr/include/pybind11 -fPIC $(python3 -m pybind11 --includes) PyODE.cpp -o _integrate$(python3-config --extension-suffix)

// #endif



















































// #ifndef PYODE_HPP
// #define PYODE_HPP


// #include <pybind11/pybind11.h>
// #include <pybind11/numpy.h>
// #include <functional>
// #include "ode.hpp"

// namespace py = pybind11;

// template<class Tt, class Ty>
// Func<Tt, Ty> to_Func(const py::object f, const std::vector<size_t>& shape);

// std::vector<size_t> shape(const py::array& arr) {
//     const ssize_t* shape_ptr = arr.shape();  // Pointer to shape data
//     size_t ndim = arr.ndim();  // Number of dimensions
//     return std::vector<size_t>(shape_ptr, shape_ptr + ndim);
// }

// template<class Tt, class Ty>
// event_f<Tt, Ty> to_event(const py::object py_event, const std::vector<size_t>& shape);

// template<class Tt, class Ty>
// is_event_f<Tt, Ty> to_event_check(const py::object py_event_check, const std::vector<size_t>& shape);

// py::dict to_PyDict(const std::map<std::string, std::vector<size_t>>& _map) {
//     py::dict py_dict;
//     for (const auto& [key, vec] : _map) {
//         py::array_t<size_t> np_array(vec.size(), vec.data()); // Create NumPy array
//         py_dict[key.c_str()] = np_array; // Assign to dictionary
//     }
//     return py_dict;
// }

// template<class Tt, class Ty>
// Ty toCPP_Array(const py::array& A);


// template<class Scalar, class ArrayType>
// py::array_t<Scalar> to_numpy(const ArrayType& array, const std::vector<size_t>& _shape = {});

// template<class Tt, class Ty>
// std::vector<Tt> flatten(const std::vector<Ty>&);

// template<class T>
// py::tuple to_tuple(const std::vector<T>& vec);

// template<class Tt, class Ty>
// std::vector<Event<Tt, Ty>> to_Events(const py::object& events, const std::vector<size_t>& shape);


// template<class Tt, class Ty>
// std::vector<StopEvent<Tt, Ty>> to_StopEvents(const py::object& events, const std::vector<size_t>& shape);


// #pragma GCC visibility push(hidden)

// template<class Tt, class Ty>
// class PyEvent : public Event<Tt, Ty>{

// public:

//     PyEvent(py::str name,
//             py::object when,
//             py::object check_if,
//             py::object mask,
//             const std::vector<size_t>& shape={})
//             :Event<Tt, Ty>(name.cast<std::string>(), to_event<Tt, Ty>(when, shape), to_event_check<Tt, Ty>(check_if, shape), to_Func<Tt, Ty>(mask, shape)), _py_when(when), _py_check_if(check_if), _py_mask(mask){}

//     PyEvent<Tt, Ty> remake(const std::vector<size_t>& shape) const {
//         return PyEvent<Tt, Ty>(this->name, _py_when, _py_check_if, _py_mask, shape);
//     }

// private:
    
//     py::object _py_when;
//     py::object _py_check_if;
//     py::object _py_mask;
// };

// template<class Tt, class Ty>
// class PyStopEvent : public StopEvent<Tt, Ty>{

// public:

//     PyStopEvent(py::str name,
//             py::object when,
//             py::object check_if,
//             const std::vector<size_t>& shape={})
//             :StopEvent<Tt, Ty>(name.cast<std::string>(), to_event<Tt, Ty>(when, shape), to_event_check<Tt, Ty>(check_if, shape)), _py_when(when), _py_check_if(check_if) {}

//     PyStopEvent<Tt, Ty> remake(const std::vector<size_t>& shape) const {
//         return PyStopEvent<Tt, Ty>(this->name, _py_when, _py_check_if, shape);
//     }

// private:
    
//     py::object _py_when;
//     py::object _py_check_if;
// };


// template<class Tt, class Ty>
// class PySolverState : public SolverState<Tt, Ty>{

// public:
//     PySolverState(const Tt& t, const Ty& q, const Tt& habs, const std::string& event, const bool& diverges, const bool& is_stiff, const bool& is_running, const bool& is_dead, const size_t& N, const std::string& message, const std::vector<size_t>& shape): SolverState<Tt, Ty>(t, q, habs, event, diverges, is_stiff, is_running, is_dead, N, message), shape(shape) {}

//     const std::vector<size_t> shape;

// };

// template<class Tt, class Ty>
// struct PyOdeResult{

//     OdeResult<Tt, Ty> src;

//     void examine() const{
//         src.examine();
//     }

// };


// template<class Tt, class Ty>
// class PyODE : public ODE<Tt, Ty>{
// public:

//     PyODE(py::object f, const Tt t0, const py::array q0, const Tt stepsize, const Tt rtol, const Tt atol, const Tt min_step, const py::tuple args, const py::str method, const Tt event_tol, py::object events, py::object stop_events) : ODE<Tt, Ty>(
//         to_Func<Tt, Ty>(f, shape(q0)), t0, toCPP_Array<Tt, Ty>(q0), stepsize, rtol, atol, min_step, 
//         toCPP_Array<Tt, std::vector<Tt>>(args), 
//         method.cast<std::string>(), event_tol, to_Events<Tt, Ty>(events, shape(q0)), to_StopEvents<Tt, Ty>(stop_events, shape(q0))), _shape(shape(q0)) {}

//     PySolverState<Tt, Ty> py_state() const{
//         SolverState<Tt, Ty> s = this->state();
//         return PySolverState<Tt, Ty>(s.t, s.q, s.habs, s.event, s.diverges, s.is_stiff, s.is_running, s.is_dead, s.N, s.message, _shape); 
//     }

//     PySolverState<Tt, Ty> py_advance(){
//         this->advance();
//         return py_state();
//     }

//     std::vector<size_t> full_shape() const{
//         std::vector<size_t> result;
//         result.reserve(1 + _shape.size()); // Pre-allocate memory for efficiency
//         result.push_back(this->t.size());        // Add the first element
//         result.insert(result.end(), _shape.begin(), _shape.end()); // Append the original vector
//         return result;
//     }

//     PyOdeResult<Tt, Ty> py_integrate(const Tt& interval, const int& max_frames, const int& max_events, const bool& terminate, const bool& display){
//         return {this->integrate(interval, max_frames, max_events, terminate, display)};
//     }

// private:
//     const std::vector<size_t> _shape;
// };



// #pragma GCC visibility pop


// /*
// -----------------------------------------------------------------------
// -----------------------------IMPLEMENTATIONS-------------------------------
// -----------------------------------------------------------------------
// */


// template<class Tt, class Ty>
// Func<Tt, Ty> to_Func(const py::object f, const std::vector<size_t>& shape){
//     if (f.is_none()){
//         return nullptr;
//     }
//     Func<Tt, Ty> g = [f, shape](const Tt& t, const Ty& y, const std::vector<Tt>& args) -> Ty {
//         return toCPP_Array<Tt, Ty>(f(t, to_numpy<Tt>(y, shape), *to_tuple(args)));
//     };
//     return g;
// }

// template<class Tt, class Ty>
// event_f<Tt, Ty> to_event(const py::object py_event, const std::vector<size_t>& shape){
//     if (py_event.is_none()){
//         return nullptr;
//     }
//     event_f<Tt, Ty> g = [py_event, shape](const Tt& t, const Ty& f, const std::vector<Tt>& args) -> Tt {
//         return py_event(t, to_numpy<Tt>(f, shape), *to_tuple(args)).template cast<Tt>();
//     };
//     return g;
// }

// template<class Tt, class Ty>
// is_event_f<Tt, Ty> to_event_check(const py::object py_event_check, const std::vector<size_t>& shape){
//     if (py_event_check.is_none()){
//         return nullptr;
//     }
//     is_event_f<Tt, Ty> g = [py_event_check, shape](const Tt& t, const Ty& f, const std::vector<Tt>& args) -> bool {
//         return py_event_check(t, to_numpy<Tt>(f, shape), *to_tuple(args)).equal(py::bool_(true));
//     };
//     return g;
// }

// template<class Tt, class Ty>
// Ty toCPP_Array(const py::array& A){
//     size_t n = A.size();
//     Ty res(n);

//     const Tt* data = static_cast<const Tt*>(A.data());

//     for (size_t i=0; i<n; i++){
//         res[i] = data[i];
//     }
//     return res;
// }

// template<class Scalar, class ArrayType>
// py::array_t<Scalar> to_numpy(const ArrayType& array, const std::vector<size_t>& _shape){
//     if (_shape.size() == 0){
//         return py::array_t<Scalar>(shape(array), array.data());
//     }
//     else{
//         return py::array_t<Scalar>(_shape, array.data());
//     }
// }

// template<class Tt, class Ty>
// std::vector<Tt> flatten(const std::vector<Ty>& f){
//     size_t nt = f.size();
//     size_t nd = f[0].size();
//     std::vector<Tt> res(nt*nd);

//     for (size_t i=0; i<nt; i++){
//         for (size_t j=0; j<nd; j++){
//             res[i*nd + j] = f[i][j];
//         }
//     }
//     return res;
// }

// template<class T>
// py::tuple to_tuple(const std::vector<T>& arr) {
//     py::tuple py_tuple(arr.size());  // Create a tuple of the same size as the vector
//     for (size_t i = 0; i < arr.size(); ++i) {
//         py_tuple[i] = py::float_(arr[i]);  // Convert each double element to py::float_
//     }
//     return py_tuple;
// }

// template<class Tt, class Ty>
// std::vector<Event<Tt, Ty>> to_Events(const py::object& events, const std::vector<size_t>& shape){   
//     if (events.is_none()){
//         return {};
//     }
//     std::vector<Event<Tt, Ty>> res;
//     for (const py::handle& item : events){
//         res.push_back(item.cast<PyEvent<Tt, Ty>&>().remake(shape));
//     }
//     return res;
// }


// template<class Tt, class Ty>
// std::vector<StopEvent<Tt, Ty>> to_StopEvents(const py::object& events, const std::vector<size_t>& shape){

//     if (events.is_none()){
//         return {};
//     }
    
//     std::vector<StopEvent<Tt, Ty>> res;
//     for (const py::handle& item : events){
//         res.push_back(item.cast<PyStopEvent<Tt, Ty>&>().remake(shape));
//     }
//     return res;
// }



// template<class Tt, class Ty>
// void define_ode_module(py::module& m) {
//     py::class_<PyEvent<Tt, Ty>>(m, "Event", py::module_local())
//         .def(py::init<py::str, py::object, py::object, py::object>(),
//             py::arg("name"),
//             py::arg("when"),
//             py::arg("check_if")=py::none(),
//             py::arg("mask")=py::none());

//     py::class_<PyStopEvent<Tt, Ty>>(m, "StopEvent", py::module_local())
//         .def(py::init<py::str, py::object, py::object>(),
//             py::arg("name"),
//             py::arg("when"),
//             py::arg("check_if")=py::none());

//     py::class_<PyOdeResult<Tt, Ty>>(m, "OdeResult", py::module_local())
//         .def("examine", &PyOdeResult<Tt, Ty>::examine);
//         // .def_property_readonly("t", [](const OdeResult<Tt, Ty>& self){
//         //     return to_numpy<Tt>(self.t);
//         // })
//         // .def_property_readonly("q", [](const OdeResult<Tt, Ty>& self){
//         //     return to_numpy<Tt>(flatten<Tt, Ty>(self.q), self.full_shape());
//         // })
//         // .def_property_readonly("events", [](const OdeResult<Tt, Ty>& self){
//         //     return to_PyDict(self.events);
//         // })
//         // .def_property_readonly("diverges", [](const OdeResult<Tt, Ty>& self){return self.diverges;})
//         // .def_property_readonly("is_stiff", [](const OdeResult<Tt, Ty>& self){return self.is_stiff;})
//         // .def_property_readonly("success", [](const OdeResult<Tt, Ty>& self){return self.success;})
//         // .def_property_readonly("runtime", [](const OdeResult<Tt, Ty>& self){return self.runtime;})
//         // .def_property_readonly("message", [](const OdeResult<Tt, Ty>& self){return self.message;})
//         // .def("examine", &OdeResult<Tt, Ty>::examine);


//     py::class_<PySolverState<Tt, Ty>>(m, "SolverState", py::module_local())
//         .def_property_readonly("t", [](const PySolverState<Tt, Ty>& self){return self.t;})
//         .def_property_readonly("q", [](const PySolverState<Tt, Ty>& self){return to_numpy<Tt>(self.q, self.shape);})
//         .def_property_readonly("event", [](const PySolverState<Tt, Ty>& self){return self.event;})
//         .def_property_readonly("diverges", [](const PySolverState<Tt, Ty>& self){return self.diverges;})
//         .def_property_readonly("is_stiff", [](const PySolverState<Tt, Ty>& self){return self.is_stiff;})
//         .def_property_readonly("is_running", [](const PySolverState<Tt, Ty>& self){return self.is_running;})
//         .def_property_readonly("is_dead", [](const PySolverState<Tt, Ty>& self){return self.is_dead;})
//         .def_property_readonly("N", [](const PySolverState<Tt, Ty>& self){return self.N;})
//         .def_property_readonly("message", [](const PySolverState<Tt, Ty>& self){return self.message;})
//         .def("show", [](const PySolverState<Tt, Ty>& self){return self.show();});

        

//     py::class_<PyODE<Tt, Ty>>(m, "LowLevelODE", py::module_local())
//         .def(py::init<py::object, Tt, py::array, Tt, Tt, Tt, Tt, py::tuple, py::str, Tt, py::object, py::object>(),
//             py::arg("f"),
//             py::arg("t0"),
//             py::arg("q0"),
//             py::arg("stepsize"),
//             py::kw_only(),
//             py::arg("rtol")=1e-6,
//             py::arg("atol")=1e-12,
//             py::arg("min_step")=0.,
//             py::arg("args")=py::tuple(),
//             py::arg("method")="RK45",
//             py::arg("event_tol")=1e-12,
//             py::arg("events")=py::none(),
//             py::arg("stop_event")=py::none())
//         .def("integrate", &PyODE<Tt, Ty>::py_integrate,
//             py::arg("interval"),
//             py::kw_only(),
//             py::arg("max_frames")=-1,
//             py::arg("max_events")=-1,
//             py::arg("terminate")=true,
//             py::arg("display")=false)
//         .def("advance", &PyODE<Tt, Ty>::py_advance)
//         .def("state", &PyODE<Tt, Ty>::py_state)
//         .def_property_readonly("t", [](const PyODE<Tt, Ty>& self){return to_numpy<Tt>(self.t);})
//         .def_property_readonly("q", [](const PyODE<Tt, Ty>& self){return to_numpy<Tt>(flatten<Tt, Ty>(self.q), self.full_shape());})
//         .def("event_map", [](const PyODE<Tt, Ty>& self){return to_PyDict(self.event_map());})
//         .def_property_readonly("runtime", [](const PyODE<Tt, Ty>& self){return self.runtime;});

// }


// //g++ -O3 -Wall -shared -std=c++20 -fopenmp -I/usr/include/python3.12 -I/usr/include/pybind11 -fPIC $(python3 -m pybind11 --includes) pyode.cpp -o _integrate$(python3-config --extension-suffix)

// #endif