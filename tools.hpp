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

template<class Tt, class Ty>
struct OdeResult{

    const std::vector<Tt> t;
    const std::vector<Ty> q;
    const std::vector<size_t> events;
    std::vector<size_t> transforms;
    bool diverges;
    bool is_stiff;
    bool success;// if t_max or break condition were satisfied
    long double runtime;
    std::string message;

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
struct SolverState{

    Tt t;
    Ty q;
    Tt habs;
    bool event;
    bool transform_event;
    bool diverges;
    bool is_stiff;
    bool is_running; //if tmax or breakcond are met or is dead, it is set to false. It can be set to true if new tmax goal is set
    bool is_dead; //e.g. if stiff or diverges. This is irreversible.
    size_t neval;
    std::string message;

    void show(const int& precision = 15) const{
        std::cout << std::endl << std::setprecision(precision) <<
        "t          : " << t << "\n" <<
        "q          : " << q.transpose() << "\n" <<
        "h          : " << habs << "\n" <<
        "event      : " << (event ? "true" : "false") << "\n" <<
        "transformed: " << (transform_event ? "true" : "false") << "\n"
        "diverges   : " << (diverges ? "true" : "false") << "\n" << 
        "stiff      : " << (is_stiff ? "true" : "false") << "\n" <<
        "running    : " << (is_running ? "true" : "false") << "\n" <<
        "dead       : " << (is_dead ? "true" : "false") << "\n" <<
        "Updates    : " << neval << "\n" <<
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


template<class Tt, class Ty, bool raw_ode, bool raw_event>
struct OdeArgs{

    SolverArgs<Tt, Ty, raw_ode, raw_event> S;
    std::string method;
};


template<class T>
std::vector<T> subvector(const std::vector<T>& a, const size_t& begin){
    size_t N = a.size();
    size_t Nsub = N - begin;

    std::vector<T> res(Nsub);
    for (size_t i=0; i<Nsub; i++){
        res[i] = a[begin+i];
    }
    return res;
}

template<class T>
std::vector<T> vec_add(const std::vector<T>& a, const T& q){
    std::vector<T> res = a;
    for (size_t i=0; i<a.size(); i++){
        res[i] += q;
    }
    return res;
}



#endif