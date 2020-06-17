from SystemDynamicsEngine import Engine
import matplotlib.pyplot as plt
import numpy as np

t = np.arange(0, 40, 0.1)
eng = Engine(t)
eng.stocks({'r': 10, 'p': 5})
eng.flow('f0', start=None, end='r', f="0.5*r")
eng.flow('f1', start='r', end=None, f="0.3*r*p")
eng.flow('f2', start=None, end='p', f="0.3*r*p")
eng.flow('f3', start='p', end=None, f="0.9*p")
eng.run()
plt.plot(t, eng.r, 'r')
plt.plot(t, eng.p, 'g')
plt.show()
