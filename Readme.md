# System Dynamics Engine in Python

This is a Python implementation of System Dynamics by Jay Forrester in Python. Additionally to the existing features coined by J. Forrester, I added new functionality with the objective of simplying the modeling process, making System Dynamics more accessible for new users in corporate environments. The objective of this simulation engine is the simulation of business models.
The System Dynamics engine in based on [JdHerman's implementation](https://github.com/jdherman/stockflow), although many changes have been made.

The backend code is running on a demo server 167.71.32.209:8000.

## Features
- Supported modeling elements: stocks, flows, clouds, variables, data sources
- Object-oriented and JSON/GUI interface for loading SD models
- Easy formula definition
- Use external data sources directly in formulas 

## Usage
Install required Python packages:
- networkx
- cexprtk
- numpy
- msgfy

Import
```python
from SystemDynamicsEngine import Engine
```
This demo will simulate a predator-prey dynamic system. First define the simulation time steps 
```python
delta_time = 0.1
t = np.arange(0, 40, delta_time)
eng = Engine(t)
```
Then define all variables used in the model
```python
eng.variables({'a': 0.5, 'b': 0.3, 'c': 0.9})
```
Next, initialize the stocks
```python
eng.stocks({'r': 10, 'p': 5})
```
The flows:
```python
eng.flow('f0', start=None, end='r', f="a*r")
eng.flow('f1', start='r', end=None, f="b*r*p")
eng.flow('f2', start=None, end='p', f="b*r*p")
eng.flow('f3', start='p', end=None, f="c*p")
```
Finally run the simulation and optionally plot the results
```python
eng.run()
plt.plot(t, eng.r, 'r')
plt.plot(t, eng.p, 'g')
plt.show()
```