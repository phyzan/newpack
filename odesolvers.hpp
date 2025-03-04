#ifndef ODESOLVERS_HPP
#define ODESOLVERS_HPP

#include <array>
#include <string>
#include "tools.hpp"
#include <limits>


template<class Tt, int N>
class OdeSolver{


public:

    using Callable = Func<Tt, N>;
    using Ty = vec<Tt, N>;

    //arguments below almost identical to solver args
    const Callable f;
    const Tt rtol;
    const Tt atol;
    const Tt h_min;
    const std::vector<Tt> args;
    const event<Tt, N> eventfunc;
    const event<Tt, N> stopevent;
    const is_event<Tt,N> check_event;
    const is_event<Tt, N> check_stop;
    const Func<Tt, N> fmask;
    const event<Tt, N> maskevent;
    const is_event<Tt, N> check_mask;
    const Tt event_tol;

    const size_t n; //size of ode system

    const Tt MAX_FACTOR = Tt(10);
    const Tt SAFETY = Tt(9)/10;
    const Tt MIN_FACTOR = Tt(2)/10;

private:

    //arguments below are passed into the SolverState when commanded
    Tt _t;
    Ty _q;
    Tt _habs;
    Tt _tmax;
    bool _event = false;
    bool _transform_event = false;
    bool _diverges = false;
    bool _is_stiff = false;
    bool _is_running = true;
    bool _is_dead = false;
    size_t _N=0;//total number of solution updates
    std::string _message = "Alive"; //different from "running".
    int _direction;


public:

    virtual ~OdeSolver() = default;

    //MODIFIERS
    void stop(const std::string& text = "") {_is_running = false; _message = (text == "") ? "Stopped by user" : text;}
    void kill(const std::string& text = "") {_is_running = false; _is_dead = true; _message = (text == "") ? "Killed by user" : text;}
    bool advance_by(const Tt& h);
    bool advance();
    void set_goal(const Tt& t_max);

    //ACCESSORS
    const Tt& t = _t;
    const Ty& q = _q;
    const Tt& stepsize = _habs;
    const Tt& tmax = _tmax;
    const int& direction = _direction;
    const bool& at_event() const {return _event;}
    const bool& at_transform_event() const {return _transform_event;}
    const bool& diverges() const {return _diverges;}
    const bool& is_stiff() const {return _is_stiff;}
    const bool& is_running() const {return _is_running;}
    const bool& is_dead() const {return _is_dead;}
    const std::string& message() {return _message;}
    const SolverState<Tt, N> state() const {
        return {_t, _q, _habs, _event, _transform_event, _diverges, _is_stiff, _is_running, _is_dead, _N, _message};
    }

    //MEMBER FUNCTIONS BELOW IMPLEMENTED BY CUSTOM DERIVED CLASSES
    //THEY MUST NOT DEPEND ON THE CURRENT STATE

    virtual Ty step(const Tt& t_old, const Ty& q_old, const Tt& h) const = 0;

    virtual State<Tt, N> adaptive_step() const = 0;



protected:

    OdeSolver(const SolverArgs<Tt, N>& S): f(S.f), rtol(S.rtol), atol(S.atol), h_min(S.h_min), args(S.args), eventfunc(S.eventfunc), stopevent(S.stopevent), check_event(S.check_event), check_stop(S.check_stop), fmask(S.fmask), maskevent(S.maskevent), check_mask(S.check_mask), event_tol(S.event_tol), n(S.q0.size()), _t(S.t0), _q(S.q0), _habs(S.habs) {
        set_goal(S.tmax);
    }


private:

    OdeSolver operator=(const OdeSolver&) = delete;
    OdeSolver(const OdeSolver& other) = default;

    bool _adapt_to_event(State<Tt, N>& next, const event<Tt, N>& event, const is_event<Tt, N>& check_event_condition)const;

    bool _go_to_state(State<Tt, N>& next);

    bool _update(const Tt& t_new, const Ty& y_new, const Tt& h_next);

    void _warn_dead(){
        throw std::runtime_error("Solver has permanently stop integrating. If this is not due to calling .kill(), call state() to see the cause.");
    }
};


/*
------------------------------------------------------------------------------
-----------------------------IMPLEMENTATIONS----------------------------------
------------------------------------------------------------------------------
*/

template<class Tt, int N>
void OdeSolver<Tt, N>::set_goal(const Tt& t_max_new){
    if ((_is_stiff || _diverges) && (!_is_dead || _is_running) ){
        //sanity check. 
        throw std::runtime_error("Bud detected");
    }

    if (_is_dead){
        _warn_dead();
    }
    else if (t_max_new == _t){
        _direction = 0;
        _tmax = t_max_new;
        stop("Direction not provided");
    }
    else{
        _tmax = t_max_new;
        _is_running = true;
        _direction = ( t_max_new > _t) ? 1 : -1;
    }
}


template<class Tt, int N>
bool OdeSolver<Tt, N>::advance(){
    State<Tt, N> next = adaptive_step();
    return _go_to_state(next);
}



template<class Tt, int N>
bool OdeSolver<Tt, N>::advance_by(const Tt& habs){
    vec<Tt, N> q_next = step(habs*direction);
    return _go_to_state({_t+habs*direction, q_next, habs});
}


template<class Tt, int N>
bool OdeSolver<Tt, N>::_update(const Tt& t_new, const vec<Tt, N>& y_new, const Tt& h_next){
    
    bool success = true;
    if (_is_dead){
        success = false;
        _warn_dead();
    }
    else if (! _is_running){
        success = false;
        throw std::runtime_error("Solver has finished integrating. Please set new t_max goal to continue integrating *before* advancing");
    }

    if (h_next < 0){//h_next is always positive
        success = false;
        throw std::runtime_error("Wrong direction of integration");
    }

    if (!All_isFinite(y_new)){
        kill();
        _diverges = true;
        _message = "Ode solution diverges";
        success = false;
    }
    else if (h_next <= h_min){
        kill();
        _is_stiff = true;
        _message = "Ode very stiff";
        success = false;
    }
    else if (t_new*direction >= _tmax*direction){
        stop();
        _q = this->step(_t, _q, _tmax-_t);
        _t = _tmax;
        _habs = h_next;
        _message = "T_max goal reached";
        _N++;
    }
    else{
        _t = t_new;
        _q = y_new;
        _habs = h_next;
        _message = "Alive";
        _N++;
    }

    return success;
}

template<class Tt, int N>
bool OdeSolver<Tt, N>::_adapt_to_event(State<Tt, N>& next, const event<Tt, N>& event, const is_event<Tt, N>& check_event_condition)const{
    // takes next state (which means tnew, hnew, and hnext_new)
    // if it is not an event or smth it is left unchanged.
    // otherwise, it is modified to depict the event with high accuracy
    if ((event != nullptr) &&  (check_event_condition == nullptr || check_event_condition(next.t, next.q, args))){
        Tt t_new, h_new;
        vec<Tt, N> q_new;

        if ( event(_t, _q, args) * event(next.t, next.q, args) <= 0 ){
            std::function<Tt(Tt)> func = [this, event](const Tt& t_next) -> Tt {
                vec<Tt, N> q_next = this->step(this->_t, this->_q, t_next-this->_t);
                return event(t_next, q_next, this->args);
            };
            
            t_new = bisect(func, _t, next.t, event_tol)[2];
            h_new = t_new - _t;
            q_new = step(_t, _q, h_new);
            next = {t_new, q_new, abs(h_new)};
            return true;
        }
        else{
            return false;
        }
    }
    else{
        return false;
    }
}

template<class Tt, int N>
bool OdeSolver<Tt, N>::_go_to_state(State<Tt, N>& next){


    if (fmask != nullptr && maskevent == nullptr){
        next.q = fmask(next.t, next.q, args);
        _transform_event = true;
    }


    if (_N > 0){
        if (_adapt_to_event(next, stopevent, check_stop)){
            bool success = _update(next.t, next.q, next.h_next);
            stop();
            _event = false;
            _transform_event = false;
            return success;
        }

        if (!_transform_event && _adapt_to_event(next, maskevent, check_mask)){
            _transform_event = true;
            _event = false;
            return _update(next.t, fmask(next.t, next.q, args), next.h_next);
        }
        else if (_adapt_to_event(next, eventfunc, check_event)){
            _event = true;
            _transform_event = false;
            return _update(next.t, next.q, next.h_next);
        }
    }
    
    _event = false;
    _transform_event = false;
    return _update(next.t, next.q, next.h_next);
}



#endif
