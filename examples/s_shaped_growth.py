from SystemDynamicsEngine import SystemDynamicsEngine
import matplotlib.pyplot as plt
import numpy as np

t = np.arange(0, 40, 0.1)
e = SystemDynamicsEngine()
e.from_json("""{ "class": "GraphLinksModel",
  "linkLabelKeysProperty": "labelKeys",
  "nodeDataArray": [ 
{"key":"p", "category":"stock", "label":"Population", "loc":"-424.15625 -44", "formula":"100", "description":"", "color":"#000000", "negative":true},
{"key":"c", "category":"stock", "label":"Customers", "loc":"-121.15625 -42", "formula":"20", "description":"", "color":"#000000", "negative":true},
{"category":"valve", "key":-3, "label":"Buying", "loc":"-272.65625 -43", "formula":"p*(0.05/600*2*0.6*c)", "description":"", "color":"#000000", "negative":true}
 ],
  "linkDataArray": [ 
{"category":"flow", "text":"flow", "from":"p", "to":"c", "labelKeys":[ -3 ]},
{"category":"influence", "text":"influence", "from":"c", "to":-3, "labelKeys":[]},
{"category":"influence", "text":"influence", "from":"p", "to":-3, "labelKeys":[]}
 ]}""")
result = e.init_sd(0, 40, 0.1, 'c', True)
plt.plot(t, result['values'])
plt.show()