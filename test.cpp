#include "ode.hpp"

using Tt = double;
using Tf = vec<Tt, 4>;

Tf f(const Tt& t, const Tf& q, const std::vector<Tt>& args){
    return {q[2], q[3], -q[0], -q[1]};
    // return q;
}

Tt fevent(const Tt& t, const Tf& f, const std::vector<Tt>& args){
    return f[1]-2;
}

Tt fevent2(const Tt& t, const Tf& f, const std::vector<Tt>& args){
    return f[1]-3;
}

bool check_fevent(const Tt& t, const Tf& f, const std::vector<Tt>& args){
    return f[3]>0;
}


int main(){

    double pi = 3.14159265359;
    double t_max = 10001*pi/2;
    Tf q0 = {1, 1, 2.3, 4.5};

    // SolverArgs<Tt, Tf, true, true> S = {f, 0., 1000, q0, 1e-3, 0., 1e-8, 0., {}, fevent, nullptr, check_fevent, nullptr, nullptr, nullptr, nullptr, 1e-12};

    Event<Tt, Tf> event1("Event1", fevent, check_fevent);
    Event<Tt, Tf> event2("Event2", fevent2, check_fevent);

    ODE<Tt, Tf> ode(f, 0, q0, 1e-2, 1e-5, 1e-10, 1e-8, {}, "RK23", 1e-10, {event1, event2});
    // ode.integrate(t_max/2, -1, 20, false);
    // ode.integrate(t_max/2, -1, 20, false);
    // OdeResultReference<Tt, Tf> res = ode.integrate(t_max, 10, 5, false);
    ode.integrate(t_max, -1, 15).examine();
    // std::cout << ode.runtime << "\n";
    // ode.state().show();
    // for (size_t i=1; i<ode.t.size(); i++){
    //     std::cout << ode.t[i]-ode.t[i-1] << "\n";
    // }
    while (true){
        ode.advance().show();
        std::cin.get();
    }

    // res.examine();

    // while (true){
        
    // }

    // res.examine();
    // ode.state().show();

    // for (size_t i=0; i<ode.t.size(); i++){
    //     std::cout << ode.t[i] << "\n";
    // }

    // for (size_t i=0; i<res.events.size(); i++){
    //     std::cout << res.q[res.events[i]].transpose() << "\n";
    // }
    // std::cout << ode.result().q[200023];

}