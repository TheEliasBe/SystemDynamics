import networkx as nx
import json
import msgfy
import cexprtk
import numpy as np
from scipy.integrate import odeint
import numpy.polynomial.polynomial as polynomial
import Intersection

class SystemDynamicsEngine:

    # graph representation
    G = nx.DiGraph()
    # table of all variables
    st = cexprtk.Symbol_Table({}, add_constants=True)
    # stores references to file handlers
    file_handlers = {}
    # stores all auxiliary variables
    aux_variables = {}
    # stores all random generators
    random_generators = {}

    # initialize custom functions
    def ite(a, b, c):
        return b if a else c
    st.functions['ite'] = ite
    def randu(a, b):
        return np.random.uniform(a,b)
    st.functions['randu'] = randu
    def randn(mu, sigma):
        return np.random.normal(mu, sigma)
    st.functions['randn'] = randn
    def randp(self, l):
        return np.random.poisson(l,1)[0]
    st.functions['randp'] = randp
    def std(*a):
        return np.std(a)


    # takes serialized graph json and returns a networkx SD graph
    def from_json(self, json_string):
        # check if input parameter is correct
        if not type(json_string) == str:
            raise Exception("Model JSON not provided as string object")
        try:
            js = json.loads(json_string)
        except Exception as exc:
            raise Exception("Could not parse JSON " + msgfy.to_error_message(exc))

        # first find all used variables and initialized them with 0
        for node in js['nodeDataArray']:
            if node['category'] == "variable" or node['category'] == 'data' or node['category'] == 'noise':
                self.st.variables[node['key']] = 0
        self.st.variables['t'] = 0
        valve_formula = {}

        # parse input
        for node in js["nodeDataArray"]:
            if node['category'] == "stock":
                if not "initialvalue" in node:
                    node["initialvalue"] = "0"
                if len(node["initialvalue"]) < 1:
                    node["initialvalue"] = "0"
                # check if formula is valid
                try:
                    self.st.variables[node['key']] = cexprtk.Expression(node["formula"], self.st).value()
                except Exception as exc:
                    raise Exception("Formula at stock '" + str(node['key']) + "' could not be parsed. " + msgfy.to_error_message(exc))
                self.G.add_node(node["key"], stock=True, formula=node["formula"], negative=node["negative"], label=node['label'], limit=100)
            elif node['category'] == "cloud":
                self.G.add_node(node["key"], cloud=True, initial_value=float("inf"))
            elif node['category'] == "valve":
                if not "formula" in node:
                    node["formula"] = "0"
                # check if formula is valid
                try:
                    cexprtk.Expression(str(node["flowformula"]), self.st).value()
                except Exception as exc:
                    raise Exception("Formula '" + str(node["flowformula"]) + "' at valve '" + str(node['key']) + "' could not be parsed. " + msgfy.to_error_message(exc))
                valve_formula[node['key']] = node['flowformula']
                self.G.add_node(node["key"], valve=True, formula=node["flowformula"], negative=False)
            elif node['category'] == "variable":
                try:
                    cexprtk.Expression(str(node["formula"]), self.st).value()
                except Exception as exc:
                    raise Exception("Formula " + str(node["formula"]) + " at variable <i>" + str(node['key']) + "</i> could not be parsed. " + msgfy.to_error_message(exc))
                self.G.add_node(node["key"], formula=node['formula'], variable=True, negative=False)
                self.st.variables[node["key"]] = float(cexprtk.Expression(node["formula"], self.st).value())
                self.aux_variables[node["key"]] = node["formula"]
            elif node['category'] == "data":
                try:
                    self.file_handlers[node['key']] =  node['src']+".csv"
                    self.st.variables[node['key']] = 0
                except Exception as ecx:
                    pass
                self.G.add_node(node["key"], data=True, source_file=node)
            elif node['category'] == "OfNodes":
                # groups have no semantic meaning when simulating
                pass
            elif node['category'] == "Comment":
                # comments have no semantic meaning when simulating
                pass
            elif node['category'].lower() == "noise":
                print("Added noise " + node['key'])
                self.random_generators[node['key']] = ("normal",0,1,node['seed'])
                self.st.variables[node['key']] = 0
                self.G.add_node(node["key"], noise=True)
                pass
            else:
                raise Exception("Undefined node type '" + node['category'] + "' in node " + str(node['label']))
        # create flows
        for flow in js["linkDataArray"]:
            if flow['category'] == 'flow':
                self.G.add_edge(flow["from"], flow["to"], formula=valve_formula[flow["labelKeys"][0]])
            elif flow['category'] == 'influence':
                # self.G.add_edge(flow["from"], flow["to"])
                pass
        return []

    # get json representation for bms tool
    def to_json(self):
        # TODO
        pass

    # create SD from networkx graph
    def from_networkx(self, networkx_model):
        self.G = networkx_model

    # returns list of all stock keys in the model as strings
    def get_stocks_str(self):
        return list(nx.get_node_attributes(self.G, 'stock').keys())

    def get_stocks_key_label_str(self):
        r = []
        for stock_key in self.get_stocks_str():
            r.append((stock_key, self.G.nodes[stock_key]["label"]))
        return r

    def get_variables_str(self):
        return list(nx.get_node_attributes(self.G, 'variable').keys())
        # return self.st.variables

    def run_simulations(self, lower_bound, upper_bound, delta_t, output_stock_keys, discrete):
        r = []
        # run simulations for each stock
        for s in output_stock_keys:
            r.append(self.init_sd(lower_bound, upper_bound, delta_t, s, discrete))
        # find intersections
        if len(output_stock_keys) > 1:
            x12 = np.arange(lower_bound, upper_bound, delta_t)
            x,y = Intersection.intersection(x12, r[0]['values'], x12, r[1]['values'])
            x = x / delta_t
            x = x + 1
            x = x.tolist()
            y = y.tolist()
            return r, x, y
        return r

    # initializes stocks and flows from the graph
    def init_sd(self, lower_bound, upper_bound, delta_t, output_stock_key, discrete):
        eng = Engine(t = np.arange(lower_bound, upper_bound, delta_t))
        # add variables to sd engine
        for variable_key in self.get_variables_str():
            print(self.G.nodes)
            print(variable_key, " ", self.G.nodes[variable_key]["formula"])
            eng.variables({variable_key: cexprtk.Expression(self.G.nodes[variable_key]["formula"], self.st).value()})

        # initialize noise generators
        for label, kind in self.random_generators.items():
            expected_value = kind[1]
            variance = kind[2]
            seed = kind[3]
            eng.noise_generators({label: (expected_value, variance, seed)})

        # add data sources to the engine variable_key:
        for label, path in self.file_handlers.items():
            eng.datasources({label: path})

        # add stocks to SD engine
        for stock_key in self.get_stocks_str():
            eng.stocks({stock_key: cexprtk.Expression(self.G.nodes[stock_key]["formula"], self.st).value()})

        # add flows
        for i,edge in enumerate(self.G.edges.data()):
            from_node = edge[0]
            to_node = edge[1]
            formula = edge[2]['formula']

            if "cloud" in self.G.nodes[from_node]:
                from_node = None
            if "cloud" in self.G.nodes[to_node]:
                to_node = None
            eng.flow(key='f'+str(i), start=from_node, end=to_node, f=formula)

        eng.run(discrete=discrete)
        run_result = list(eng.get_result(output_stock_key))
        x_max_index = np.where(run_result == np.amax(run_result))[0]
        y_max = float(np.amax(run_result))
        x_min_index = np.where(run_result == np.amin(run_result))[0].tolist()
        y_min = float(np.amin(run_result))
        x_max_index = x_max_index + 1
        x_max_index = x_max_index.tolist()
        return {"label": self.G.nodes[output_stock_key]["label"], "values": run_result, "y_max": y_max, "x_max_index": x_max_index, "y_min": y_min, "x_min_index": x_min_index}


class Engine:


    def __init__(self, t):
        # simulation time steps as np-array
        self.t = t
        # reference to all items to make sure keys are unique
        self.ix = {}
        # reference to each flow
        self.flows = {}
        # reference to each stock
        self.stock = {}
        # stocks to be simulated
        self.current = []
        # flag to see if simulation was run
        self.done = False
        # stores the simulated value of each stock
        self.results = None
        # table of all variables
        self.st = cexprtk.Symbol_Table({}, add_constants=True)
        # set to start time
        self.st.variables['t'] = 0
        # stores all data source variables (variable_key, formula)
        self.data = {}
        # stores all random number sequences for noise generators
        self.noise = {}

        def random_series(t, s):
            np.random.seed(int(t) + int(s))
            return np.random.normal()

        self.st.functions['random'] = random_series


    def flow(self, key, f, start=None, end=None):
        # resolve any data sources in formula
        for k, v in self.data.items():
            if k in f:
                f = f.replace(k, v)

        for k, v  in self.noise.items():
            if k in f:
                f = f.replace(k, "random(t,"+str(v)+")")

        # resolve noise in formula

        self.__new_state_var(key, cexprtk.Expression(f, self.st).value())
        s = self.ix[start] if start is not None else None
        e = self.ix[end] if end is not None else None
        # add flows to dict
        self.flows[key] = {'f': f, 'start': s, 'end': e}

    def f(self, y, t):
        self.current = y # set y to the initial condition
        d = np.zeros((len(y),))
        for k, f in self.flows.items():  # calculate flows only once. distribute to stocks.
            i = self.ix[k]
            # update time
            self.st.variables['t'] = t
            # update all stock variables
            for k,v in self.stock.items():
                self.st.variables[k] = self.get_stock_value(k)
            # update noise generators
            current_time_step_index = np.where(self.t == t)[0]


            ft = cexprtk.Expression(f['f'], self.st).value()
            d[i] = ft - self.current[i]
            if f['start'] is not None: d[f['start']] -= ft
            if f['end'] is not None: d[f['end']] += ft
        return d

    def run(self, discrete=False):
        self.done = False
        if not discrete:
            # odeint(f(y,t), initial condition, sequence of points)
            self.results = odeint(self.f, self.current, self.t)
        else:
            self.results = np.zeros((len(self.t), len(self.current)))
            self.results[0, :] = self.current
            for i in np.arange(1, len(self.t)):
                self.results[i, :] = self.results[i - 1, :] + self.f(self.results[i - 1, :], self.t[i])
        self.done = True
        self.current = self.results[0, :]  # restore initial conditions

    def stocks(self, icdict):
        for k, v in icdict.items():
            # initialize stock values as variables in formulas
            if type(v) == int or type(v) == float:
                self.st.variables[k] = v
            elif type(v) == str:
                self.st.variables[k] = cexprtk.Expression(v, self.st).value()
            else:
                pass
            self.__new_state_var(k, v)
            self.stock[k] = v

    def variables(self, icdict):
        self.st.variables['rand'] = 0
        for k, v in icdict.items():
            # initializes variables as variables in formulas
            self.st.variables[k] = v

    def datasources(self, icdict):
        # key (label) value (source file)
        for k, v in icdict.items():
            data_np_array = np.genfromtxt(v)
            self.data[k] = self.data_interpolation(data_np_array)

    def noise_generators(self, icdict):
        # key (label) value (source file)
        for k, v in icdict.items():
            expected_value = v[0]
            variance = v[1]
            seed = v[2]
            self.noise[k] = seed

    def __getattr__(self, key):
        if not self.done:
            return self.current[self.ix[key]]
        else:
            return self.results[:, self.ix[key]]

    def get_result(self, key):
        if not self.done:
            return self.current[self.ix[key]]
        else:
            return self.results[:, self.ix[key]]

    def get_stock_value(self, key):
        return self.current[self.ix[key]]

    def __validate_key(self, key):
        if key in self.ix: raise NameError("Variable " + key + " already defined.")

    def __new_state_var(self, key, IC):
        self.__validate_key(key)
        self.current.append(IC)
        self.ix[key] = len(self.current) - 1

    def data_interpolatio(self, source_array, size, method):
        result_array = np.zeros(size)
        if method == "duplicate":
            for i in range(size):
                data_index = int(np.linspace(0, len(source_array)-1, size)[i])
                result_array[i] = source_array[data_index]
        elif method == "polyint":
            x = np.linspace(0, len(source_array), len(source_array))
            x_new = np.linspace(0, len(source_array), size)
            data_poly = polynomial.polyfit(x=x, y=source_array, deg=max(3,len(source_array)//2))
            result_array = polynomial.polyval(x_new, data_poly)
        elif method == "padlast":
            for i in range(size):
                result_array[i] = source_array[-1]
            pass
        elif method == "padmean":
            mean = np.mean(source_array)
            for i in range(size):
                if i < len(source_array):
                    result_array[i] = source_array[i]
                else:
                    result_array[i] = mean
            pass
        elif method == "padzero":
            for i in range(len(source_array)):
                result_array[i] = source_array[i]
        elif method == "reflect":
            # TODO
            pass
        x = np.linspace(0, len(result_array)-1, len(result_array))
        result_coeff = polynomial.polyfit(x=x, y=result_array, deg=10)
        result_coeff_repr = []
        for i, c in enumerate(result_coeff):
            result_coeff_repr.append(str(c)+"*t^"+str(i))
        return "+".join(result_coeff_repr)

    def data_interpolation(self, source_array):
        # stretch interpolation nodes to new interval if necessary
        x = np.linspace(0, self.t[-1], len(source_array))
        # apply polynomial interpolation
        data_poly = polynomial.polyfit(x=x, y=source_array, deg=min(10,max(3, len(source_array) // 2)))
        # create formula as string
        data_poly = np.round(data_poly, 4)
        result_coeff_repr = []
        for i, c in enumerate(data_poly):
            result_coeff_repr.append(str(c)+"*t^"+str(i))
        interpolation_formula = '+'.join(result_coeff_repr)
        return interpolation_formula





