from SystemDynamicsEngine import Engine
import matplotlib.pyplot as plt
import numpy as np

t = np.arange(0, 40, 0.1)
eng = Engine(t)
eng.stocks({'r': 10, 'p': 5})
eng.datasources({'data': 'example.csv'})
eng.flow('f0', start=None, end='r', f="data")
eng.run()
plt.plot(eng.r, 'r')
plt.show()
