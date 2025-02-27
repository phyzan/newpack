#ifndef ODE_HPP
#define ODE_HPP

#include <variant>
#include "rk_adaptive.hpp"
#include <unordered_map>
#include <chrono>
#include <omp.h>


template<class Tt, class Ty, bool raw_ode, bool raw_event>
OdeSolver<Tt, Ty, raw_ode, raw_event>* getSolver(const OdeArgs<Tt, Ty, raw_ode, raw_event>& args) {

    OdeSolver<Tt, Ty, raw_ode, raw_event>* solver = nullptr;

    if (args.method == "RK23") {
        solver = new RK23<Tt, Ty, raw_ode, raw_event>(args.S);
    }
    else if (args.method == "RK45") {
        solver = new RK45<Tt, Ty, raw_ode, raw_event>(args.S);
    }
    else {
        throw std::runtime_error("Unknown solver method");
    }

    return solver;
}

template<class Tt, class Ty, bool raw_ode = true, bool raw_event = true>
class ODE{

public:

    ODE(const OdeArgs<Tt, Ty, raw_ode, raw_event>& params) : _solver(getSolver(params)) {
        _t_arr.push_back(params.S.t0);
        _q_arr.push_back(params.S.q0);
    }

    ~ODE(){delete _solver;}

    const OdeResult<Tt, Ty> integrate(const Tt& interval, const int& max_frames=-1, const int& max_events=-1, const bool& terminate = true, const bool& display = false);

    const SolverState<Tt, Ty> state() const {return _solver->state();}

    const std::vector<size_t>& events() const{return _events;}
    const std::vector<size_t>& transforms() const{return _transforms;}
    const long double& runtime() {return _runtime;}
    const std::vector<Ty>& q() const{ return _q_arr;}
    const std::vector<Tt>& t() const{return _t_arr;}



private:

    OdeSolver<Tt, Ty, raw_ode, raw_event>* _solver;
    std::vector<Tt> _t_arr;
    std::vector<Ty> _q_arr;
    std::vector<size_t> _events;
    std::vector<size_t> _transforms;
    long double _runtime = 0.;
};


/*
-----------------------------------------------------------------------
-----------------------------IMPLEMENTATIONS-------------------------------
-----------------------------------------------------------------------
*/


template<class Tt, class Ty, bool raw_ode, bool raw_event>
const OdeResult<Tt, Ty> ODE<Tt, Ty, raw_ode, raw_event>::integrate(const Tt& interval, const int& max_frames, const int& max_events, const bool& terminate, const bool& display){

    if (interval <= 0){
        throw std::runtime_error("Enter a positive interval");
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    const Tt t0 = _solver->t_now();
    const size_t N = _t_arr.size();
    const size_t Nevents = _events.size();
    const size_t Ntransfomrs = _transforms.size();
    long int event_counter = 0;
    long int frame_counter = 0;
    size_t i = N;
    Tt t;
    Ty q;

    _solver->set_goal(t0+interval*_solver->direction);

    while (_solver->is_running()){
        if (_solver->advance()){
            t = _solver->t_now();
            q = _solver->q_now();
            if (_solver->at_event() && event_counter < max_events){
                _events.push_back(i);
                _t_arr.push_back(t);
                _q_arr.push_back(q);
                if (++event_counter == max_events && terminate){
                    _solver->stop("Max events reached");
                }
                ++i;
            }
            else if ( (max_frames == -1) || (abs(t-t0)*max_frames > frame_counter*interval) ){
                _t_arr.push_back(t);
                _q_arr.push_back(q);
                if (_solver->at_transform_event() && _solver->maskevent != nullptr){
                    _transforms.push_back(i);
                }
                ++frame_counter;
                // if (++frame_counter == max_frames) _solver->stop("Max frames reached"); //shouldnt be required
                ++i;
            }
        }
    }

    OdeResult<Tt, Ty> res = {subvector(_t_arr, N), subvector(_q_arr, N), vec_add(subvector(_events, Nevents), -N), vec_add(subvector(_transforms, Ntransfomrs), -N), _solver->diverges(), _solver->is_stiff(), !_solver->is_dead(), 0, _solver->message()};

    auto t2 = std::chrono::high_resolution_clock::now();
    

    std::chrono::duration<long double> rt = t2-t1;

    res.runtime = rt.count();
    _runtime += res.runtime;
    return res;
}


#endif