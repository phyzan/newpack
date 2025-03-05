from _integrate import *
from scipy.integrate import solve_ivp
import time
import numpy as np


def f(t, q):
    return [q[2], q[3], -q[0], -q[1]]

def event(t, q):
    return q[1]

def check_event(t, q):
    return q[3] > 0


ps_event = Event("Poincare Section", when=event, check_if=check_event, mask=f)

w = StopEvent("Poincare Section", when=event, check_if=check_event)


q0 = np.array([1, 1, 2.3, 4.5])


ode = LowLevelODE(f, 0., q0, 0.01, atol=1e-8, rtol=0., event_tol=0., events=[ps_event], stop_events=[w])

res = ode.integrate(1000).examine()
print(ode.runtime)
# ode.state().show()
# while True:
#     ode.state().show()
#     input()
#     ode.advance()