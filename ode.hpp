#ifndef ODE_HPP
#define ODE_HPP

#include <variant>
#include <span>
#include "rk_adaptive.hpp"
#include "rk_classic.hpp"
#include <unordered_map>
#include <chrono>
#include <omp.h>


template<class Tt, class Ty>
OdeSolver<Tt, Ty>* getSolver(const SolverArgs<Tt, Ty>& S, const std::string& method) {

    OdeSolver<Tt, Ty>* solver = nullptr;

    if (method == "RK23") {
        solver = new RK23<Tt, Ty>(S);
    }
    else if (method == "RK45") {
        solver = new RK45<Tt, Ty>(S);
    }
    else if (method == "RK4") {
        solver = new RK4<Tt, Ty>(S);
    }
    else {
        throw std::runtime_error("Unknown solver method");
    }

    return solver;
}

template<class Tt, class Ty>
class ODE{

public:

    ODE(const Func<Tt, Ty> f, const Tt t0, const Ty q0, const Tt stepsize, const Tt rtol, const Tt atol, const Tt min_step, const std::vector<Tt> args = {}, const std::string& method = "RK45", const Tt event_tol = 1e-10, const std::vector<Event<Tt, Ty>>& events = {}, const std::vector<StopEvent<Tt, Ty>>& stop_events = {}) : _event_holder(events) {

        const SolverArgs<Tt, Ty> S = {f, t0, t0, q0, stepsize, rtol, atol, min_step, args, events, stop_events, event_tol};
        _solver = getSolver(S, method);
        _register_state();
    }

    ODE(const SolverArgs<Tt, Ty>& S, const std::string& method) : _event_holder(S.events){
        _solver = getSolver(S, method);
        _register_state();
    }

    ~ODE(){delete _solver;}

    const OdeResultReference<Tt, Ty> integrate(const Tt& interval, const int& max_frames=-1, const int& max_events=-1, const bool& terminate = true, const bool& display = false);

    const SolverState<Tt, Ty> state() const {return _solver->state();}

    SolverState<Tt, Ty> advance();

    const std::vector<Tt>& t = _t_arr;
    const std::vector<Ty>& q = _q_arr;
    const long double& runtime = _runtime;

private:

    EventHolder<Tt, Ty> _event_holder;
    OdeSolver<Tt, Ty>* _solver;
    std::vector<Tt> _t_arr;
    std::vector<Ty> _q_arr;
    long double _runtime = 0.;

    void _register_state(){
        _t_arr.push_back(_solver->t);
        _q_arr.push_back(_solver->q);
    }
};






/*
-----------------------------------------------------------------------
-----------------------------IMPLEMENTATIONS-------------------------------
-----------------------------------------------------------------------
*/


template<class Tt, class Ty>
const OdeResultReference<Tt, Ty> ODE<Tt, Ty>::integrate(const Tt& interval, const int& max_frames, const int& max_events, const bool& terminate, const bool& display){
    
    auto t1 = std::chrono::high_resolution_clock::now();

    const Tt t0 = _solver->t;
    const size_t N = _t_arr.size();
    long int event_counter = 0;
    long int frame_counter = 0;
    size_t i = N;

    _solver->set_goal(t0+interval);
    if (max_events == 0){
        _solver->stop("Max events was set to 0. No integration performed");
    }


    while (_solver->is_running()){
        if (_solver->advance()){

            if ((event_counter != max_events) && _solver->at_event()){
                _event_holder.include(*(_solver->current_event()), _t_arr.size());
                if ( (++event_counter == max_events) && terminate){
                    _solver->stop("Max events reached");
                }
                ++i;
                _register_state();
            }
            else if ( (max_frames == -1) || (abs(_solver->t-t0)*max_frames >= (frame_counter+1)*interval) ){
                _register_state();
                ++frame_counter;
                ++i;
            }
        }
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<long double> rt = t2-t1;

    OdeResultReference<Tt, Ty> res = {std::span<Tt>(_t_arr).subspan(N), std::span<Ty>(_q_arr).subspan(N), _event_holder.subHolder(N), _solver->diverges(), _solver->is_stiff(), !_solver->is_dead(), rt.count(), _solver->message()};

    _runtime += res.runtime;
    return res;
}


template<class Tt, class Ty>
SolverState<Tt, Ty> ODE<Tt, Ty>::advance(){
    if (!_solver->is_running()){
        _solver->set_goal(std::numeric_limits<Tt>::infinity());
    }
    if (_solver->advance()){
        if (_solver->current_event() != nullptr){
            _event_holder.include(*_solver->current_event(), _t_arr.size());
        }
        _register_state();
    }
    
    _solver->set_goal(_solver->t);
    return _solver->state();
}


#endif