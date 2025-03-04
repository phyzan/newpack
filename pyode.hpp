#ifndef PYODE_HPP
#define PYODE_HPP


#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <functional>
#include "ode.hpp"

namespace py = pybind11;

template<class Tt, class Ty>
Func<Tt, Ty> to_Func(const py::object f, const std::vector<size_t>& shape);

std::vector<size_t> shape(const py::array& arr) {
    const ssize_t* shape_ptr = arr.shape();  // Pointer to shape data
    size_t ndim = arr.ndim();  // Number of dimensions
    return std::vector<size_t>(shape_ptr, shape_ptr + ndim);
}

template<class Tt, class Ty>
event_f<Tt, Ty> to_event(const py::object py_event, const std::vector<size_t>& shape);

template<class Tt, class Ty>
is_event_f<Tt, Ty> to_event_check(const py::object py_event_check, const std::vector<size_t>& shape);

template<class Tt, class Ty>
Ty toCPP_Array(const py::array& A);


template<class Scalar, class ArrayType>
py::array_t<Scalar> to_numpy(const ArrayType& array, const std::vector<size_t>& _shape = {});

template<class Tt, class Ty>
std::vector<Tt> flatten(const std::vector<Ty>&);

template<class T>
py::tuple to_tuple(const std::vector<T>& vec);


#pragma GCC visibility push(hidden)
template<class Tt, class Ty>
struct PyOdeResult{
    
    const OdeResult<Tt, Ty> res_src;

    py::array_t<Tt> t() const {
        return to_numpy<Tt>(res_src.t);
    }

    py::array_t<Tt> q() const {
        return to_numpy<Tt>(flatten<Tt, Ty>(res_src.q), res_src.full_shape());
    }

    py::array_t<size_t> events() const {
        return to_numpy<size_t>(res_src.events);
    }

    py::array_t<size_t> transforms() const {
        return to_numpy<size_t>(res_src.transforms);
    }

    bool diverges() const{
        return res_src.diverges;
    }

    bool is_stiff() const{
        return res_src.is_stiff;
    }

    bool success() const{
        return res_src.success;
    }

    long double runtime() const{
        return res_src.runtime;
    }

    py::str message() const{
        return py::str(res_src.message);
    }

    void examine() const{
        res_src.examine();
    }

};


template<class Tt, class Ty>
struct PySolverState{
    const SolverState<Tt, Ty> ss;

    const Tt t() const {
        return ss.t;
    }

    const py::array_t<Tt> q() const {
        return to_numpy<Tt>(ss.q);
    }

    bool event() const {
        return ss.event;
    }

    bool transform_event() const {
        return ss.transform_event;
    }

    bool diverges() const{
        return ss.diverges;
    }

    bool is_stiff() const{
        return ss.is_stiff;
    }

    bool is_running() const{
        return ss.is_running;
    }

    bool is_dead() const{
        return ss.is_dead;
    }

    size_t N() const{
        return size_t(ss.N);
    }

    py::str message() const{
        return py::str(ss.message);
    }

    void show() const{
        ss.show();
    }
};


template<class Tt, class Ty>
struct PyOdeArgs{

    const py::object f;
    const Tt t0;
    const py::array q0;
    const Tt stepsize;
    const Tt rtol;
    const Tt atol;
    const Tt min_step;
    const py::tuple args;
    const py::str method;
    const Tt event_tol;
    const py::object event;
    const py::object stopevent;
    const py::object check_event;
    const py::object check_stop;
    const py::object fmask;
    const py::object maskevent;
    const py::object check_mask;

    ODE<Tt, Ty, false, false> to_ODE() const {
        std::vector<size_t> s = shape(q0);
        return ODE<Tt, Ty, false, false>(to_Func<Tt, Ty>(f, s), t0, toCPP_Array<Tt, Ty>(q0), stepsize, rtol, atol, min_step, toCPP_Array<Tt, std::vector<Tt>>(args), method.cast<std::string>(), event_tol, to_event<Tt, Ty>(event, s), to_event<Tt, Ty>(stopevent, s), to_event_check<Tt, Ty>(check_event, s), to_event_check<Tt, Ty>(check_stop, s), to_Func<Tt, Ty>(fmask, s), to_event<Tt, Ty>(maskevent, s), to_event_check<Tt, Ty>(check_mask, s));
    }

};


template<class Tt, class Ty>
class PyODE{

private:

    ODE<Tt, Ty, false, false> ode;
    std::vector<size_t> _shape;

    std::vector<size_t> fullshape() {
        std::vector<size_t> result;
        result.reserve(1 + _shape.size()); // Pre-allocate memory for efficiency
        result.push_back(ode.t.size());        // Add the first element
        result.insert(result.end(), _shape.begin(), _shape.end()); // Append the original vector
        return result;
    }

public:

    PyODE(const py::object f, const Tt& t0, const py::array& q0, const Tt& stepsize, const Tt& rtol, const Tt& atol, const Tt& min_step, const py::tuple& args, const py::str& method, const Tt& event_tol, const py::object event, const py::object stopevent, const py::object check_event, const py::object check_stop, const py::object fmask, const py::object maskevent, const py::object check_mask) : ode(PyOdeArgs<Tt, Ty>{f, t0, q0, stepsize, rtol, atol, min_step, args, method, event_tol, event, stopevent, check_event, check_stop, fmask, maskevent, check_mask}.to_ODE()), _shape(shape(q0)) {}

    PyODE(const PyOdeArgs<Tt, Ty>& S) : ode(S.to_ODE()) {}

    const PyOdeResult<Tt, Ty> integrate(const Tt& interval, const int& max_frames=-1, const int& max_events=-1, const bool& terminate = true, const bool& display = false){
        return {ode.integrate(interval, max_frames, max_events, terminate, display).as_copy()};
    }

    const PySolverState<Tt, Ty> state() const {
        return {ode.state()};
    }

    PySolverState<Tt, Ty> advance(){
        return {ode.advance()};
    }

    const py::array_t<Tt> t(){return to_numpy<Tt>(ode.t);}
    const py::array_t<Tt> q(){return to_numpy<Tt>(flatten<Tt, Ty>(ode.q), fullshape());}
    const py::array_t<size_t> events(){return to_numpy<size_t>(ode.events);}
    const py::array_t<size_t> transforms(){return to_numpy<size_t>(ode.transforms);}
    const long double runtime(){return ode.runtime;}

};


#pragma GCC visibility pop


/*
-----------------------------------------------------------------------
-----------------------------IMPLEMENTATIONS-------------------------------
-----------------------------------------------------------------------
*/


template<class Tt, class Ty>
Func<Tt, Ty> to_Func(const py::object f, const std::vector<size_t>& shape){
    if (f.is_none()){
        return nullptr;
    }
    Func<Tt, Ty> g = [f, shape](const Tt& t, const Ty& y, const std::vector<Tt>& args) -> Ty {
        return toCPP_Array<Tt, Ty>(f(t, to_numpy<Tt>(y, shape), *to_tuple(args)));
    };
    return g;
}

template<class Tt, class Ty>
event_f<Tt, Ty> to_event(const py::object py_event, const std::vector<size_t>& shape){
    if (py_event.is_none()){
        return nullptr;
    }
    event_f<Tt, Ty> g = [py_event, shape](const Tt& t, const Ty& f, const std::vector<Tt>& args) -> Tt {
        return py_event(t, to_numpy<Tt>(f, shape), *to_tuple(args)).template cast<Tt>();
    };
    return g;
}

template<class Tt, class Ty>
is_event_f<Tt, Ty> to_event_check(const py::object py_event_check, const std::vector<size_t>& shape){
    if (py_event_check.is_none()){
        return nullptr;
    }
    is_event_f<Tt, Ty> g = [py_event_check, shape](const Tt& t, const Ty& f, const std::vector<Tt>& args) -> bool {
        return py_event_check(t, to_numpy<Tt>(f, shape), *to_tuple(args)).equal(py::bool_(true));
    };
    return g;
}

template<class Tt, class Ty>
Ty toCPP_Array(const py::array& A){
    size_t n = A.size();
    Ty res(n);

    const Tt* data = static_cast<const Tt*>(A.data());

    for (size_t i=0; i<n; i++){
        res[i] = data[i];
    }
    return res;
}

template<class Scalar, class ArrayType>
py::array_t<Scalar> to_numpy(const ArrayType& array, const std::vector<size_t>& _shape){
    if (_shape.size() == 0){
        return py::array_t<Scalar>(shape(array), array.data());
    }
    else{
        return py::array_t<Scalar>(_shape, array.data());
    }
}

template<class Tt, class Ty>
std::vector<Tt> flatten(const std::vector<Ty>& f){
    size_t nt = f.size();
    size_t nd = f[0].size();
    std::vector<Tt> res(nt*nd);

    for (size_t i=0; i<nt; i++){
        for (size_t j=0; j<nd; j++){
            res[i*nd + j] = f[i][j];
        }
    }
    return res;
}

template<class T>
py::tuple to_tuple(const std::vector<T>& arr) {
    py::tuple py_tuple(arr.size());  // Create a tuple of the same size as the vector
    for (size_t i = 0; i < arr.size(); ++i) {
        py_tuple[i] = py::float_(arr[i]);  // Convert each double element to py::float_
    }
    return py_tuple;
}



template<class Tt, class Ty>
void define_ode_module(py::module& m) {
    py::class_<PyODE<Tt, Ty>>(m, "LowLevelODE", py::module_local())
        .def(py::init<py::object, Tt, py::array, Tt, Tt, Tt, Tt, py::tuple, py::str, Tt, py::object, py::object, py::object, py::object, py::object, py::object, py::object>(),
            py::arg("f"),
            py::arg("t0"),
            py::arg("q0"),
            py::arg("stepsize"),
            py::kw_only(),
            py::arg("rtol")=1e-6,
            py::arg("atol")=1e-12,
            py::arg("min_step")=0.,
            py::arg("args")=py::tuple(),
            py::arg("method")="RK45",
            py::arg("event_tol")=1e-12,
            py::arg("event")=py::none(),
            py::arg("stopevent")=py::none(),
            py::arg("check_event")=py::none(),
            py::arg("check_stop")=py::none(),
            py::arg("fmask")=py::none(),
            py::arg("maskevent")=py::none(),
            py::arg("check_mask")=py::none())
        .def("integrate", &PyODE<Tt, Ty>::integrate,
            py::arg("interval"),
            py::kw_only(),
            py::arg("max_frames")=-1,
            py::arg("max_events")=-1,
            py::arg("terminate")=true,
            py::arg("display")=false)
        .def("advance", &PyODE<Tt, Ty>::advance)
        .def("state", &PyODE<Tt, Ty>::state)
        .def_property_readonly("t", &PyODE<Tt, Ty>::t)
        .def_property_readonly("q", &PyODE<Tt, Ty>::q)
        .def_property_readonly("events", &PyODE<Tt, Ty>::events)
        .def_property_readonly("transforms", &PyODE<Tt, Ty>::transforms)
        .def_property_readonly("runtime", &PyODE<Tt, Ty>::runtime);


    py::class_<PyOdeResult<Tt, Ty>>(m, "OdeResult", py::module_local())
        // .def(py::init<py::array, py::array, bool, bool, double>(), py::arg("t"), py::arg("y"), py::arg("diverges"), py::arg("is_stiff"), py::arg("runtime"){})
        .def_property_readonly("t", &PyOdeResult<Tt, Ty>::t)
        .def_property_readonly("q", &PyOdeResult<Tt, Ty>::q)
        .def_property_readonly("events", &PyOdeResult<Tt, Ty>::events)
        .def_property_readonly("transforms", &PyOdeResult<Tt, Ty>::transforms)
        .def_property_readonly("diverges", &PyOdeResult<Tt, Ty>::diverges)
        .def_property_readonly("is_stiff", &PyOdeResult<Tt, Ty>::is_stiff)
        .def_property_readonly("success", &PyOdeResult<Tt, Ty>::success)
        .def_property_readonly("runtime", &PyOdeResult<Tt, Ty>::runtime)
        .def_property_readonly("message", &PyOdeResult<Tt, Ty>::message)
        .def("examine", &PyOdeResult<Tt, Ty>::examine);

    py::class_<PySolverState<Tt, Ty>>(m, "SolverState", py::module_local())
    .def_property_readonly("t", &PySolverState<Tt, Ty>::t)
    .def_property_readonly("q", &PySolverState<Tt, Ty>::q)
    .def_property_readonly("is_event", &PySolverState<Tt, Ty>::event)
    .def_property_readonly("is_transform_event", &PySolverState<Tt, Ty>::transform_event)
    .def_property_readonly("diverges", &PySolverState<Tt, Ty>::diverges)
    .def_property_readonly("is_stiff", &PySolverState<Tt, Ty>::is_stiff)
    .def_property_readonly("is_running", &PySolverState<Tt, Ty>::is_running)
    .def_property_readonly("is_dead", &PySolverState<Tt, Ty>::is_dead)
    .def_property_readonly("N", &PySolverState<Tt, Ty>::N)
    .def_property_readonly("message", &PySolverState<Tt, Ty>::message)
    .def("show", &PySolverState<Tt, Ty>::show);
}


//g++ -O3 -Wall -shared -std=c++20 -fopenmp -I/usr/include/python3.12 -I/usr/include/pybind11 -fPIC $(python3 -m pybind11 --includes) PyODE.cpp -o _integrate$(python3-config --extension-suffix)

#endif