from _integrate import *
from scipy.integrate import solve_ivp
import time
import numpy as np

def f(t, q):
    return [q[2], q[3], -q[0], -q[1]]

def event(t, q):
    return q[1]


ode = LowLevelODE(f, 0., np.array([1, 1, 2.3, 4.5]), 0.01, event_tol=0., event=event)

ode.integrate(10001)
ode.integrate(-10001)
print(ode.q[-1])
# while True:
#     ode.state().show()
#     input()
#     ode.advance()