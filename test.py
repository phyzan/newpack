from _integrate import *
from scipy.integrate import solve_ivp
import time
import numpy as np

def f(t, q):
    return np.array([[q[1, 0], q[1, 1]], [-q[0, 0], -q[0, 1]]])
    return [q[2], q[3], -q[0], -q[1]]

def event(t, q):
    return q[0, 1]

q0 = np.array([1, 1, 2.3, 4.5])
q0 = np.array([[1, 1], [2.3, 4.5]])


ode = LowLevelODE(f, 0., q0, 0.01, event_tol=0., event=event)

ode.integrate(1001)
ode.integrate(-1001)
ode.q[:, :] = 0
print(ode.q is ode.q)
# while True:
#     ode.state().show()
#     input()
#     ode.advance()