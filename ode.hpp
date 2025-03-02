#ifndef ODE_HPP
#define ODE_HPP

#include <variant>
#include <span>
#include "rk_adaptive.hpp"
#include "rk_classic.hpp"
#include <unordered_map>
#include <chrono>
#include <omp.h>


template<class Tt, class Ty, bool raw_ode, bool raw_event>
OdeSolver<Tt, Ty, raw_ode, raw_event>* getSolver(const SolverArgs<Tt, Ty, raw_ode, raw_event>& S, const std::string& method) {

    OdeSolver<Tt, Ty, raw_ode, raw_event>* solver = nullptr;

    if (method == "RK23") {
        solver = new RK23<Tt, Ty, raw_ode, raw_event>(S);
    }
    else if (method == "RK45") {
        solver = new RK45<Tt, Ty, raw_ode, raw_event>(S);
    }
    else if (method == "RK4") {
        solver = new RK4<Tt, Ty, raw_ode, raw_event>(S);
    }
    else {
        throw std::runtime_error("Unknown solver method");
    }

    return solver;
}

template<class Tt, class Ty, bool raw_ode = true, bool raw_event = true>
class ODE{

public:

    ODE(const ode_t<Tt, Ty, raw_ode> f, const Tt t0, const Ty q0, const Tt stepsize, const Tt rtol, const Tt atol, const Tt min_step, const std::vector<Tt> args = {},  const std::string& method = "RK45", const Tt event_tol = 1e-10, const event_t<Tt, Ty, raw_event> event = nullptr, const event_t<Tt, Ty, raw_event> stopevent = nullptr, const is_event_t<Tt, Ty, raw_event> check_event = nullptr, const is_event_t<Tt, Ty, raw_event> check_stop = nullptr, const mask<Tt, Ty> fmask = nullptr, const event_f<Tt, Ty> maskevent = nullptr, const is_event_f<Tt, Ty> check_mask = nullptr) {

        const SolverArgs<Tt, Ty, raw_ode, raw_event> S = {f, t0, t0, q0, stepsize, rtol, atol, min_step, args, event, stopevent, check_event, check_stop, fmask, maskevent, check_mask, event_tol};
        _solver = getSolver(S, method);
        _t_arr.push_back(S.t0);
        _q_arr.push_back(S.q0);
    }

    ODE(const SolverArgs<Tt, Ty, raw_ode, raw_event>& S, const std::string& method){
        _solver = getSolver(S, method);
        _t_arr.push_back(S.t0);
        _q_arr.push_back(S.q0);
    }

    ~ODE(){delete _solver;}

    const OdeResultReference<Tt, Ty> integrate(const Tt& interval, const int& max_frames=-1, const int& max_events=-1, const bool& terminate = true, const bool& display = false);

    const SolverState<Tt, Ty> state() const {return _solver->state();}

    const std::vector<size_t> _recent(const std::vector<size_t>& event_list, const size_t& begin, const size_t& Tbegin) const{

        const size_t N = event_list.size()-begin;
        std::vector<size_t> res(N);
        for (size_t i=0; i<N; i++){
            res[i] = event_list[begin+i]-Tbegin;
        }
        return res;

    }

    SolverState<Tt, Ty> advance();

    const std::vector<Tt>& t = _t_arr;
    const std::vector<Ty>& q = _q_arr;
    const std::vector<size_t>& events = _events;
    const std::vector<size_t>& transforms = _transforms;
    const long double& runtime = _runtime;

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
const OdeResultReference<Tt, Ty> ODE<Tt, Ty, raw_ode, raw_event>::integrate(const Tt& interval, const int& max_frames, const int& max_events, const bool& terminate, const bool& display){
    
    auto t1 = std::chrono::high_resolution_clock::now();

    const Tt t0 = _solver->t;
    const size_t N = _t_arr.size();
    const size_t Nevents = _events.size();
    const size_t Ntransfomrs = _transforms.size();
    long int event_counter = 0;
    long int frame_counter = 0;
    size_t i = N;
    Tt _t;
    Ty _q;

    _solver->set_goal(t0+interval);

    while (_solver->is_running()){
        if (_solver->advance()){
            _t = _solver->t;
            _q = _solver->q;
            if (_solver->at_event() && event_counter != max_events){
                _events.push_back(i);
                _t_arr.push_back(_t);
                _q_arr.push_back(_q);
                if (++event_counter == max_events && terminate){
                    _solver->stop("Max events reached");
                }
                ++i;
            }
            else if ( (max_frames == -1) || (abs(_t-t0)*max_frames >= (frame_counter+1)*interval) ){
                _t_arr.push_back(_t);
                _q_arr.push_back(_q);
                if (_solver->at_transform_event() && _solver->maskevent != nullptr){
                    _transforms.push_back(i);
                }
                ++frame_counter;
                ++i;
            }
        }
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<long double> rt = t2-t1;

    OdeResultReference<Tt, Ty> res = {std::span<Tt>(_t_arr).subspan(N), std::span<Ty>(_q_arr).subspan(N), _recent(_events, Nevents, N), _recent(_transforms, Ntransfomrs, N), _solver->diverges(), _solver->is_stiff(), !_solver->is_dead(), rt.count(), _solver->message()};

    _runtime += res.runtime;
    return res;
}


template<class Tt, class Ty, bool raw_ode, bool raw_event>
SolverState<Tt, Ty> ODE<Tt, Ty, raw_ode, raw_event>::advance(){
    if (!_solver->is_running()){
        _solver->set_goal(std::numeric_limits<Tt>::infinity());
    }
    _solver->advance();
    _t_arr.push_back(_solver->t);
    _q_arr.push_back(_solver->q);
    size_t i = _t_arr.size()-1;
    if (_solver->at_event()){
        _events.push_back(i);
    }
    else if (_solver->at_transform_event() && _solver->maskevent != nullptr){
        _transforms.push_back(i);
    }
    
    _solver->set_goal(_solver->t);
    return _solver->state();
}


#endif