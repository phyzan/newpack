#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <vector>
#include <iostream>
#include <iomanip>
#include <eigen3/Eigen/Dense>


// USEFUL ALIASES

template<class Tt, size_t Nr=0, size_t Nc=1>
using vec = std::conditional_t<(Nr == 0), Eigen::Array<Tt, Eigen::Dynamic, Nc>, Eigen::Array<Tt, Nr, Nc>>;

template<class Tt, class Ty>
using ode = Ty(*)(const Tt&, const Ty&, const std::vector<Tt>&);

template<class Tt, class Ty>
using ode_f = std::function<Ty(const Tt&, const Ty&, const std::vector<Tt>&)>;

template<class Tt, class Ty, bool raw>
using ode_t = std::conditional_t<(raw==true), ode<Tt, Ty>, ode_f<Tt, Ty>>;



template<class Tt, class Ty>
using event = Tt(*)(const Tt&, const Ty&, const std::vector<Tt>&);

template<class Tt, class Ty>
using event_f = std::function<Tt(const Tt&, const Ty&, const std::vector<Tt>&)>;

template<class Tt, class Ty, bool raw>
using event_t = std::conditional_t<(raw==true), event<Tt, Ty>, event_f<Tt, Ty>>;


template<class Tt, class Ty>
using is_event = bool(*)(const Tt&, const Ty&, const std::vector<Tt>&);

template<class Tt, class Ty>
using is_event_f = std::function<bool(const Tt&, const Ty&, const std::vector<Tt>&)>;

template<class Tt, class Ty, bool raw>
using is_event_t = std::conditional_t<(raw==true), is_event<Tt, Ty>, is_event_f<Tt, Ty>>;



template<class Tt, class Ty>
using mask = ode_f<Tt, Ty>; //mask function is of the same form as the ode: f(t, y) -> y_new


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
std::vector<int> shape(const Eigen::Array<T, Nr, Nc>& arr){
    return {arr.rows(), arr.cols()};
}

template<class T>
std::vector<int> shape(const std::vector<T>& arr){
    return {arr.size()};
}

//BISECTION USED FOR EVENTS IN ODES

template<class T, typename Callable>
std::vector<T> bisect(Callable&& f, const T& a, const T& b, const T& atol){
    T err = 2*atol;
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


//ODERESULT STRUCT TO ENCAPSULATE THE RESULT OF AN ODE INTEGRATION

template<class Tt, class Ty, template <class> class Container>
class _OdeResult{

public:

    const Container<Tt> t;
    const Container<Ty> q;
    const std::vector<size_t> events;
    const std::vector<size_t> transforms;
    const bool diverges;
    const bool is_stiff;
    const bool success;// if t_max or break condition were satisfied
    const long double runtime;
    const std::string message;

    void examine() const{
        std::cout << std::endl <<
        "Points           : " << t.size() << "\n" <<
        "Events           : " << events.size() << "\n" <<
        "Transformations  : " << transforms.size() << "\n"
        "Diverges         : " << (diverges ? "true" : "false") << "\n" << 
        "Stiff            : " << (is_stiff ? "true" : "false") << "\n" <<
        "Success          : " << (success ? "true" : "false") << "\n" <<
        "Runtime          : " << runtime << "\n" <<
        "Termination cause: " << message << "\n";
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
                this->events, this->transforms, this->diverges, this->is_stiff, this->success, this->runtime, this->message};
    }
};


template<class Tt, class Ty>
struct SolverState{

    const Tt t;
    const Ty q;
    const Tt habs;
    const bool event;
    const bool transform_event;
    const bool diverges;
    const bool is_stiff;
    const bool is_running; //if tmax or breakcond are met or is dead, it is set to false. It can be set to true if new tmax goal is set
    const bool is_dead; //e.g. if stiff or diverges. This is irreversible.
    const size_t N;
    const std::string message;

    void show(const int& precision = 15) const{
        std::cout << std::endl << std::setprecision(precision) <<
        "t          : " << t << "\n" <<
        "q          : " << q.transpose() << "\n" <<
        "h          : " << habs << "\n" <<
        "Event      : " << (event ? "true" : "false") << "\n" <<
        "Transformed: " << (transform_event ? "true" : "false") << "\n"
        "Diverges   : " << (diverges ? "true" : "false") << "\n" << 
        "Stiff      : " << (is_stiff ? "true" : "false") << "\n" <<
        "Running    : " << (is_running ? "true" : "false") << "\n" <<
        "Updates    : " << N << "\n" <<
        "Dead       : " << (is_dead ? "true" : "false") << "\n" <<
        "State      : " << message << "\n";
    }

};


template<class Tt, class Ty>
struct State{

    Tt t;
    Ty q;
    Tt h_next;
};


template<class Tt, class Ty, bool raw_ode, bool raw_event>
struct SolverArgs{

    const ode_t<Tt, Ty, raw_ode> f;
    const Tt t0;
    const Tt tmax;
    const Ty q0;
    const Tt habs;
    const Tt rtol;
    const Tt atol;
    const Tt h_min;
    const std::vector<Tt> args;
    const event_t<Tt, Ty, raw_event> event;
    const event_t<Tt, Ty, raw_event> stopevent;
    const is_event_t<Tt, Ty, raw_event> check_event;
    const is_event_t<Tt, Ty, raw_event> check_stop;
    const mask<Tt, Ty> fmask;
    const event_f<Tt, Ty> maskevent;
    const is_event_f<Tt, Ty> check_mask;
    const Tt event_tol;
};


#endif