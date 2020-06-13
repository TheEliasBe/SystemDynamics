# System Dynamics Engine in Python

This is a Python implementation of System Dynamics by Jay Forrester in Python. Additionally to the existing features coined by J. Forrester, I added new functionality with the objective of simplying the modeling process, making System Dynamics more accessible for new users in corporate environments. The objective of this simulation engine is the simulation of business models.

The backend code is running on a demo server 167.71.32.209:8000.

## Usage
Install required Python packages:
- networkx
- cexprtk
- numpy
- matplotlib
- msgfy

Clone repository into your project

Import your model
```python
from SystemDynamicsEngine import SystemDynamicsEngine
e = SystemDynamicsEngine()
e.from_json(my_model_as_json_string)
e.eval_stock("stock1", lowerbound=0, upperbound=10)
```

