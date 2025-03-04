#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <vector>
#include <iostream>
#include <iomanip>
#include <eigen3/Eigen/Dense>


// USEFUL ALIASES

template<class Tt, int N>
using vec = Eigen::Array<Tt, 1, N>;

template<class Tt, int N>
using VecArray = Eigen::Array<Tt, -1, N>;

template<class Tt, int N>
using Func = std::function<vec<Tt, N>(const Tt&, const vec<Tt, N>&, const std::vector<Tt>&)>;

template<class Tt, int N>
using event = std::function<Tt(const Tt&, const vec<Tt, N>&, const std::vector<Tt>&)>;

template<class Tt, int N>
using is_event = std::function<bool(const Tt&, const vec<Tt, N>&, const std::vector<Tt>&)>;


template<class T, int Nr, int Nc>
T norm(const Eigen::Array<T, Nr, Nc>& f){
    return (f*f).sum();
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


//ODERESULT STRUCT TO ENCAPSULATE THE RESULT OF AN ODE INTEGRATION

template<class tContainer, class yContainer>
class _OdeResult{

public:

    const tContainer t;
    const yContainer q;
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

template<class Tt, int N>
using OdeResult = _OdeResult<std::vector<Tt>, VecArray<Tt, N>>;

template<class Tt, int N>
class OdeResultReference : public _OdeResult<std::span<Tt>, Eigen::Map<VecArray<Tt, N>>> {

public:
    OdeResult<Tt, N> as_copy() const{
        return {std::vector<Tt>(this->t.begin(), this->t.end()),
                this->q.eval(),
                this->events, this->transforms, this->diverges, this->is_stiff, this->success, this->runtime, this->message};
    }
};


template<class Tt, int N>
struct SolverState{
    
    const Tt t;
    const vec<Tt, N> q;
    const Tt habs;
    const bool event;
    const bool transform_event;
    const bool diverges;
    const bool is_stiff;
    const bool is_running; //if tmax or breakcond are met or is dead, it is set to false. It can be set to true if new tmax goal is set
    const bool is_dead; //e.g. if stiff or diverges. This is irreversible.
    const size_t Nt;
    const std::string message;

    void show(const int& precision = 15) const{
        std::cout << std::endl << std::setprecision(precision) <<
        "t          : " << t << "\n" <<
        "q          : " << q << "\n" <<
        "h          : " << habs << "\n" <<
        "Event      : " << (event ? "true" : "false") << "\n" <<
        "Transformed: " << (transform_event ? "true" : "false") << "\n"
        "Diverges   : " << (diverges ? "true" : "false") << "\n" << 
        "Stiff      : " << (is_stiff ? "true" : "false") << "\n" <<
        "Running    : " << (is_running ? "true" : "false") << "\n" <<
        "Updates    : " << Nt << "\n" <<
        "Dead       : " << (is_dead ? "true" : "false") << "\n" <<
        "State      : " << message << "\n";
    }

};


template<class Tt, int N>
struct State{

    Tt t;
    vec<Tt, N> q;
    Tt h_next;
};


template<class Tt, int N>
struct SolverArgs{

    const Func<Tt, N> f;
    const Tt t0;
    const Tt tmax;
    const vec<Tt, N> q0;
    const Tt habs;
    const Tt rtol;
    const Tt atol;
    const Tt h_min;
    const std::vector<Tt> args;
    const event<Tt, N> eventfunc;
    const event<Tt, N> stopevent;
    const is_event<Tt, N> check_event;
    const is_event<Tt, N> check_stop;
    const Func<Tt, N> fmask;
    const event<Tt, N> maskevent;
    const is_event<Tt, N> check_mask;
    const Tt event_tol;
};


#endif