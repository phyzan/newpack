#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <vector>
#include <iostream>
#include <map>
#include <iomanip>
#include <eigen3/Eigen/Dense>


// USEFUL ALIASES

template<class Tt, int N>
using vec = Eigen::Array<Tt, 1, N>;

template<class Tt, class Ty>
using Func = std::function<Ty(const Tt&, const Ty&, const std::vector<Tt>&)>;

template<class Tt, class Ty>
using event_f = std::function<Tt(const Tt&, const Ty&, const std::vector<Tt>&)>;

template<class Tt, class Ty>
using is_event_f = std::function<bool(const Tt&, const Ty&, const std::vector<Tt>&)>;



template<class T, int Nr, int Nc>
T norm(const Eigen::Array<T, Nr, Nc>& f){
    return (f*f).sum();
}

template<class T, int Nr, int Nc>
Eigen::Array<T, Nr, Nc> cwise_abs(const Eigen::Array<T, Nr, Nc>& f){
    return f.cwiseAbs();
}

template<class T, int Nr, int Nc>
Eigen::Array<T, Nr, Nc> cwise_max(const Eigen::Array<T, Nr, Nc>& a, const Eigen::Array<T, Nr, Nc>& b){
    return a.cwiseMax(b);
}

template<class T>
T abs(const T& x){
    return (x > 0) ? x : -x;
}

template<class T, int Nr, int Nc>
bool All_isFinite(const Eigen::Array<T, Nr, Nc>& arr){
    return arr.isFinite().all();
}

template<class T, int Nr, int Nc>
std::vector<size_t> shape(const Eigen::Array<T, Nr, Nc>& arr){
    return {arr.rows(), arr.cols()};
}

template<class T>
std::vector<size_t> shape(const std::vector<T>& arr){
    return {arr.size()};
}

//BISECTION USED FOR EVENTS IN ODES

template<class T, typename Callable>
std::vector<T> bisect(Callable&& f, const T& a, const T& b, const T& atol){
    T err = 2*atol+1;
    T _a = a;
    T _b = b;
    T c = a;
    T fm;

    if (f(a)*f(b) > 0){
        throw std::runtime_error("Root not bracketed");
    }

    while (err > atol){
        c = (_a+_b)/2;
        if (c == _a || c == _b){
            break;
        }
        fm = f(c);
        if (f(_a) * fm  > 0){
            _a = c;
        }
        else{
            _b = c;
        }
        err = abs(fm);
    }
    return {_a, c, _b};
}



template<class Tt, class Ty>
class Event{

public:
    Event(const std::string& name, event_f<Tt, Ty> when, is_event_f<Tt, Ty> check_if=nullptr, Func<Tt, Ty> mask=nullptr) : _name(name), _when(when), _check_if(check_if), _mask(mask) {}

    Event(const std::string& name, std::vector<Tt> when, is_event_f<Tt, Ty> check_if=nullptr, Func<Tt, Ty> mask=nullptr) : _name(name), _when_t(when), _check_if(check_if), _mask(mask) {}

    bool determine(const Tt& t1, const Tt& t2, const std::vector<Tt>& args, std::function<Ty(const Tt&)> q, const Tt& tol){
        delete _t_event;
        delete _q_event;
        _t_event = nullptr;
        _q_event = nullptr;
        if (_when != nullptr){
            std::function<Tt(Tt)> obj_fun = [this, q, args](const Tt& t) ->Tt {
                return _when(t, q(t), args);
            };

            if (_check_if == nullptr || (_check_if(t1, q(t1), args) && _check_if(t2, q(t2), args))){
                if (_when(t1, q(t1), args) * _when(t2, q(t2), args) <= 0){

                    _t_event = new Tt;
                    _q_event = new Ty;
                    *_t_event = bisect(obj_fun, t1, t2, tol)[2];
                    *_q_event = (_mask == nullptr) ? q(*_t_event) : _mask(*_t_event, q(*_t_event), args);
                    return true;
                }
            }
        }
        else if ( (t1 <= _when_t[_i]) && (_when_t[_i] <= t2)){
            _t_event = new Tt;
            _q_event = new Ty;
            *_t_event = _when_t[_i++];
            *_q_event = (_mask == nullptr) ? q(*_t_event) : _mask(*_t_event, q(*_t_event), args);
            return true;
        }
        return false;
    }

    const Tt& t_event(){
        return *_t_event;
    }

    const Ty& q_event(){
        return *_q_event;
    }

    const std::string& name = _name;

    ~Event(){
        delete _t_event;
        delete _q_event;
    }


private:

    const std::string _name;
    const std::vector<Tt> _when_t = {};
    const event_f<Tt, Ty> _when;
    const is_event_f<Tt, Ty> _check_if = nullptr;
    const Func<Tt, Ty> _mask = nullptr;
    int _i = 0;
    Tt* _t_event = nullptr;
    Ty* _q_event = nullptr;
};





template<class Tt, class Ty>
class StopEvent{

public:
    StopEvent(const std::string& name, event_f<Tt, Ty> when, is_event_f<Tt, Ty> check_if) : _name(name), _when(when), _check_if(check_if) {}

    bool is_between(const Tt& t1, const Ty& q1, const Tt& t2, const Ty& q2, const std::vector<Tt>& args)const{
        if (_check_if == nullptr || (_check_if(t1, q1, args) && _check_if(t2, q2, args))){
            return _when(t1, q1, args) * _when(t2, q2, args) <= 0;
        }
        else{
            return false;
        }
    }

    const std::string& name = _name;


private:

    const std::string _name;
    const event_f<Tt, Ty> _when = nullptr;
    const is_event_f<Tt, Ty> _check_if = nullptr;

};


//ODERESULT STRUCT TO ENCAPSULATE THE RESULT OF AN ODE INTEGRATION

template<class Tt, class Ty, template <class> class Container>
class _OdeResult{

public:

    const Container<Tt> t;
    const Container<Ty> q;
    const std::map<std::string, std::vector<size_t>> events;
    const bool diverges;
    const bool is_stiff;
    const bool success;// if the OdeSolver didnt die during the integration
    const long double runtime;
    const std::string message;

    void examine() const{
        std::cout << std::endl <<
        "Points           : " << t.size();
        examine_events(); std::cout <<
        "Diverges         : " << (diverges ? "true" : "false") << "\n" << 
        "Stiff            : " << (is_stiff ? "true" : "false") << "\n" <<
        "Success          : " << (success ? "true" : "false") << "\n" <<
        "Runtime          : " << runtime << "\n" <<
        "Termination cause: " << message << "\n";
    }

    std::vector<size_t> full_shape() const{
        std::vector<size_t> result;
        std::vector<size_t> _shape = shape(q[0]);
        result.reserve(1 + _shape.size()); // Pre-allocate memory for efficiency
        result.push_back(t.size());        // Add the first element
        result.insert(result.end(), _shape.begin(), _shape.end()); // Append the original vector
        return result;
    }

    void examine_events() const{
        std::cout << std::endl <<
        "Events:\n----------\n";
        for (const auto& [name, array] : events){
            std::cout << "    " << name << " : " << array.size() << "\n";
        }
        std::cout << "\n----------\n";
    }
    
};

template<class Tt, class Ty>
using OdeResult = _OdeResult<Tt, Ty, std::vector>;

template<class Tt, class Ty>
class OdeResultReference : public _OdeResult<Tt, Ty, std::span> {

public:
    OdeResult<Tt, Ty> as_copy() const{
        return {std::vector<Tt>(this->t.begin(), this->t.end()),
                std::vector<Ty>(this->q.begin(), this->q.end()),
                this->events, this->diverges, this->is_stiff, this->success, this->runtime, this->message};
    }
};


template<class Tt, class Ty>
struct SolverState{

    const Tt t;
    const Ty q;
    const Tt habs;
    const std::map<std::string, bool> is_event;
    const bool diverges;
    const bool is_stiff;
    const bool is_running; //if tmax or breakcond are met or is dead, it is set to false. It can be set to true if new tmax goal is set
    const bool is_dead; //e.g. if stiff or diverges. This is irreversible.
    const size_t N;
    const std::string message;

    void show(const int& precision = 15) const{
        std::cout << std::endl << std::setprecision(precision) <<
        "t          : " << t << "\n" <<
        "q          : " << q << "\n" <<
        "h          : " << habs << "\n";
        show_event_state(); std::cout <<
        "Diverges   : " << (diverges ? "true" : "false") << "\n" << 
        "Stiff      : " << (is_stiff ? "true" : "false") << "\n" <<
        "Running    : " << (is_running ? "true" : "false") << "\n" <<
        "Updates    : " << N << "\n" <<
        "Dead       : " << (is_dead ? "true" : "false") << "\n" <<
        "State      : " << message << "\n";
    }


    void show_event_state() const{
        std::cout << std::endl <<
        "Events:\n----------\n";
        for (const auto& [name, value] : is_event) {
            std::cout << "    " << name << ": " << ( value ? "true" : "false") << "\n";
        }
        std::cout << "\n----------\n";
    }

};


template<class Tt, class Ty>
struct State{

    Tt t;
    Ty q;
    Tt h_next;
};


template<class Tt, class Ty>
struct SolverArgs{

    const Func<Tt, Ty> f;
    const Tt t0;
    const Tt tmax;
    const Ty q0;
    const Tt habs;
    const Tt rtol;
    const Tt atol;
    const Tt h_min;
    const std::vector<Tt> args;
    const std::vector<Event<Tt, Ty>> events;
    const std::vector<StopEvent<Tt, Ty>> stop_events;
    const Tt event_tol;
};









#endif