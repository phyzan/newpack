#ifndef ODE_HPP
#define ODE_HPP

#include <variant>
#include <span>
#include "rk_adaptive.hpp"
#include <unordered_map>
#include <chrono>
#include <omp.h>


template<class Tt, int N>
OdeSolver<Tt, N>* getSolver(const SolverArgs<Tt, N>& S, const std::string& method) {

    OdeSolver<Tt, N>* solver = nullptr;

    if (method == "RK23") {
        solver = new RK23<Tt, N>(S);
    }
    else if (method == "RK45") {
        solver = new RK45<Tt, N>(S);
    }
    // else if (method == "RK4") {
    //     solver = new RK4<Tt, N>(S);
    // }
    else {
        throw std::runtime_error("Unknown solver method");
    }

    return solver;
}

template<class Tt, int N>
class ODE{


public:


    ODE(const Func<Tt, N> f, const Tt t0, const vec<Tt, N> q0, const Tt stepsize, const Tt rtol, const Tt atol, const Tt min_step, const std::vector<Tt> args = {},  const std::string& method = "RK45", const Tt event_tol = 1e-10, const event<Tt, N> eventfunc = nullptr, const event<Tt, N> stopevent = nullptr, const is_event<Tt, N> check_event = nullptr, const is_event<Tt, N> check_stop = nullptr, const Func<Tt, N> fmask = nullptr, const event<Tt, N> maskevent = nullptr, const is_event<Tt, N> check_mask = nullptr) {

        const SolverArgs<Tt, N> S = {f, t0, t0, q0, stepsize, rtol, atol, min_step, args, eventfunc, stopevent, check_event, check_stop, fmask, maskevent, check_mask, event_tol};
        _solver = getSolver(S, method);
        _append(t0, q0);
    }

    ODE(const SolverArgs<Tt, N>& S, const std::string& method){
        _solver = getSolver(S, method);
        _append(S.t0, S.q0);
    }

    ~ODE(){delete _solver;}

    const OdeResult<Tt, N> integrate(const Tt& interval, const int& max_frames=-1, const int& max_events=-1, const bool& terminate = true, const bool& display = false);

    const SolverState<Tt, N> state() const {return _solver->state();}

    const std::vector<size_t> _recent(const std::vector<size_t>& event_list, const size_t& begin, const size_t& Tbegin) const{

        const size_t Nt = event_list.size()-begin;
        std::vector<size_t> res(Nt);
        for (size_t i=0; i<Nt; i++){
            res[i] = event_list[begin+i]-Tbegin;
        }
        return res;

    }

    SolverState<Tt, N> advance();

    const std::vector<Tt>& t = _t_arr;
    const VecArray<Tt, N>& q = _q_arr;
    const std::vector<size_t>& events = _events;
    const std::vector<size_t>& transforms = _transforms;
    const long double& runtime = _runtime;

private:

    OdeSolver<Tt, N>* _solver;
    std::vector<Tt> _t_arr;
    VecArray<Tt, N> _q_arr;
    std::vector<size_t> _events;
    std::vector<size_t> _transforms;
    long double _runtime = 0.;

    void _append(const Tt& tnew, const vec<Tt, N>& qnew){
        _t_arr.push_back(tnew);
        _q_arr.conservativeResize(_q_arr.rows()+1, Eigen::NoChange);
        _q_arr.row(_q_arr.rows() - 1) = qnew;
    }
};






/*
-----------------------------------------------------------------------
-----------------------------IMPLEMENTATIONS-------------------------------
-----------------------------------------------------------------------
*/


template<class Tt, int N>
const OdeResult<Tt, N> ODE<Tt, N>::integrate(const Tt& interval, const int& max_frames, const int& max_events, const bool& terminate, const bool& display){
    
    auto t1 = std::chrono::high_resolution_clock::now();

    const Tt t0 = _solver->t;
    const size_t Nt = _t_arr.size();
    const size_t Nevents = _events.size();
    const size_t Ntransfomrs = _transforms.size();
    long int event_counter = 0;
    long int frame_counter = 0;
    size_t i = Nt;
    Tt _t;
    vec<Tt, N> _q;

    _solver->set_goal(t0+interval);

    while (_solver->is_running()){
        if (_solver->advance()){
            _t = _solver->t;
            _q = _solver->q;
            if (_solver->at_event() && event_counter != max_events){
                _events.push_back(i);
                _append(_t, _q);
                if (++event_counter == max_events && terminate){
                    _solver->stop("Max events reached");
                }
                ++i;
            }
            else if ( (max_frames == -1) || (abs(_t-t0)*max_frames >= (frame_counter+1)*interval) ){
                _append(_t, _q);
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
    // OdeResult<Tt, N> res = {std::vector<Tt>(_t_arr.data()+Nt, _t_arr.data()+_t_arr.size()), _q_arr.matrix().transpose().block(Nt, 0, _q_arr.rows()-Nt, _solver->n).array(), _recent(_events, Nevents, Nt), _recent(_transforms, Ntransfomrs, Nt), _solver->diverges(), _solver->is_stiff(), !_solver->is_dead(), rt.count(), _solver->message()};

    OdeResult<Tt, N> res = {_t_arr, _q_arr, _events, _transforms, _solver->diverges(), _solver->is_stiff(), !_solver->is_dead(), rt.count(), _solver->message()};

    _runtime += res.runtime;
    return res;
}


template<class Tt, int N>
SolverState<Tt, N> ODE<Tt, N>::advance(){
    if (!_solver->is_running()){
        _solver->set_goal(std::numeric_limits<Tt>::infinity());
    }
    _solver->advance();
    _append(_solver->t, _solver->q);
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