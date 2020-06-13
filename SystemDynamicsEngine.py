import networkx as nx
import json
import msgfy
import cexprtk
import numpy as np
from scipy.integrate import odeint
import numpy.polynomial.polynomial as polynomial

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


    # takes serialized graph json from from GoJS and returns a SD graph
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
            if node['category'] == "variable" or node['category'] == 'data':
                self.st.variables[node['key']] = 0
        self.st.variables['t'] = 0
        valve_formula = {}

        # parse input
        for node in js["nodeDataArray"]:
            if node['category'] == "stock":
                if not  "formula" in node:
                    node["formula"] = "0"
                if not  "negative" in node:
                    node["negative"] = True
                if len(node["formula"]) < 1:
                    node["formula"] = "0"
                # check if formula is valid
                try:
                    self.st.variables[node['key']] = cexprtk.Expression(node["formula"], self.st).value()
                except Exception as exc:
                    raise Exception("Formula at stock '" + str(node['key']) + "' could not be parsed. " + msgfy.to_error_message(exc))
                self.G.add_node(node["key"], stock=True, formula=node["formula"], negative=node["negative"], label=node['label'], limit=100)
            elif node['category'] == "cloud":
                self.G.add_node(node["key"], initial_value=float("inf"))
            elif node['category'] == "valve":
                if not "formula" in node:
                    node["formula"] = "0"
                # check if formula is valid
                try:
                    cexprtk.Expression(str(node["formula"]), self.st).value()
                except Exception as exc:
                    raise Exception("Formula '" + str(node["formula"]) + "' at valve '" + str(node['key']) + "' could not be parsed. " + msgfy.to_error_message(exc))
                valve_formula[node['key']] = node['formula']
                self.G.add_node(node["key"], valve=True, formula=node["formula"], negative=False)
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
                    self.file_handlers[node['key']] = np.genfromtxt("user_provided_data/" + node['src']+".csv")
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
            elif node['category'] == "noise":
                # TODO
                pass
            else:
                raise Exception("Undefined node type '" + node['category'] + "' in node " + str(node['label']))
        # create flows
        for flow in js["linkDataArray"]:
            if flow['category'] == 'flow':
                self.G.add_edge(flow["from"], flow["to"], formula=valve_formula[flow["labelKeys"][0]])
            elif flow['category'] == 'influence':
                self.G.add_edge(flow["from"], flow["to"])
        return []


    def to_json(self):
        # TODO
        pass


    def from_networkx(self, networkx_model):
        self.G = networkx_model


    def data_resizer(self, source_array, size, method):
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
        print(len(x), len(result_array))
        result_coeff = polynomial.polyfit(x=x, y=result_array, deg=len(result_array)-1)
        result_coeff_repr = []
        for i, c in enumerate(result_coeff):
            result_coeff_repr.append(str(c)+"*t^"+str(i))
        return "+".join(result_coeff_repr)


    def compute_flow_analytically(self, stock_key, t_0, upper_bound, delta_t):
        predecessors = list(self.G.predecessors(stock_key))
        successors = list(self.G.successors(stock_key))
        result_values = []
        number_of_increments = int(upper_bound / delta_t)

        # set initial value for stock
        result_values.append(cexprtk.Expression(self.G.nodes[stock_key]["formula"], self.st).value())

        # re-size data sources
        resized_data_arrays = {}
        for key in self.file_handlers.keys():
            resized_data_arrays[key] = self.data_resizer(self.file_handlers[key], number_of_increments, "polyint")

        # resolve auxiliary variables
        for s in successors + predecessors:
            for v in self.aux_variables.keys():
                if v in self.G.nodes[s]["formula"]:
                    self.G.nodes[s]["formula"] = self.G.nodes[s]["formula"].replace(v, self.G.nodes[v]['formula'])

        # set noise generators
        for s in successors + predecessors:
            for v in self.random_generators.keys():
                self.st.variables[v] = np.random.normal(1,0)

        for i in range(number_of_increments)[1:]:
            self.st.variables['t'] = t_0 + i * delta_t
            inflow_value = 0
            outflow_value = 0
            # update external data values
            for d in self.file_handlers.keys():
                self.st.variables[d] = resized_data_arrays[d][i]

            # set value of stock as variables
            self.st.variables[stock_key] = result_values[i-1]
            for k in successors:
                k = next(self.G.successors(k))
                self.st.variables[k] = self.compute_flow_analytically(k, t_0, t_0 + i * delta_t, delta_t)[1][-1]

            # accumulated inflow
            for p in predecessors:
                predecessors_inflow = cexprtk.Expression(self.G.nodes[p]["formula"], self.st).value()
                predecessor_stock_key = next(self.G.predecessors(p))
                predecessor_stock_value = self.compute_flow_analytically(predecessor_stock_key, t_0, t_0 + i * delta_t, delta_t)

                # preceding stock is allowed to go negative, just add the inflow
                if self.G.nodes[predecessor_stock_key]["negative"] == True:
                    inflow_value += predecessors_inflow
                else:
                # preceeding stock is not allowed to go negative
                # check if stock would pass 0
                    if predecessor_stock_value[0] == True:
                        # preceding stock went negative, no inflow from that stock
                        inflow_value += 0
                    elif predecessor_stock_value[1][-1] >= predecessors_inflow:
                        # preceding stock is positive and value greater than flow
                        inflow_value += predecessors_inflow
                    elif predecessor_stock_value[1][-1] < predecessors_inflow:
                        # preceding stock is positive but value is smaller than flow -> limit flow to current value
                        inflow_value += predecessor_stock_value[1][-1]

            # accumulate outflow
            for s in successors:
                outflow_value += cexprtk.Expression(self.G.nodes[s]["formula"], self.st).value()

            # compute the net inflow to this stock
            net_flow = result_values[i-1] + inflow_value*delta_t - outflow_value*delta_t

            # check if this stock is allowed to be negative
            if self.G.nodes[stock_key]["negative"]:
                result_values.append(net_flow)
            else:
                if net_flow > 0:
                    result_values.append(net_flow)
                else:
                    return (True, result_values)

        return (False, result_values)

    # returns list of all stock keys in the model as strings
    def get_stocks_str(self):
        return list(nx.get_node_attributes(self.G, 'stock').keys())

    def function_evaluator(self, t, expr_str):
        self.st.variables['t'] = t
        return cexprtk.Expression(expr_str, self.st).value()

    # initializes stocks and flows from the graph
    def init_sd(self, lower_bound, upper_bound, delta_t, output_stock_key):
        eng = Engine(t = np.arange(lower_bound, upper_bound, delta_t))
        # add stocks to SD engine
        for stock_key in self.get_stocks_str():
            eng.stocks({stock_key: cexprtk.Expression(self.G.nodes[stock_key]["formula"], self.st).value()})

        # add flows
        for i,edge in enumerate(self.G.edges.data()):
            from_node = edge[0]
            to_node = edge[1]
            formula = edge[2]['formula']
            eng.flow(key='f'+str(i), start=from_node, end=to_node, f=formula)

        eng.run()
        return [{"label": output_stock_key, "values": list(eng.get_result(output_stock_key))}]


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


    def flow(self, key, f, start=None, end=None):
        self.__new_state_var(key, cexprtk.Expression(f, self.st).value())
        s = self.ix[start] if start is not None else None
        e = self.ix[end] if end is not None else None
        # resolve any data sources in formula
        for k, v in self.data.items():
            if k in f:
                f = f.replace(k, v)
        # add flows to dict
        self.flows[key] = {'f': f, 'start': s, 'end': e}

    def xdot(self, y, t):
        self.current = y
        d = np.zeros((len(y),))
        for k, f in self.flows.items():  # calculate flows only once. distribute to stocks.
            i = self.ix[k]
            # update time
            self.st.variables['t'] = t
            # update all stock variables
            for k,v in self.stock.items():
                self.st.variables[k] = self.get_stock_value(k)
            ft = cexprtk.Expression(f['f'], self.st).value()
            d[i] = ft - self.current[i]
            if f['start'] is not None: d[f['start']] -= ft
            if f['end'] is not None: d[f['end']] += ft
        return d

    def run(self, discrete=False):
        self.done = False
        if not discrete:
            self.results = odeint(self.xdot, self.current, self.t)
        else:
            self.results = np.zeros((len(self.t), len(self.current)))
            self.results[0, :] = self.current
            for i in np.arange(1, len(self.t)):
                self.results[i, :] = self.results[i - 1, :] + self.xdot(self.results[i - 1, :], self.t[i])
        self.done = True
        self.current = self.results[0, :]  # restore initial conditions

    def stocks(self, icdict):
        for k, v in icdict.items():
            # initialize stock values as variables in formulas
            self.st.variables[k] = v
            self.__new_state_var(k, v)
            self.stock[k] = v

    def variables(self, icdict):
        for k, v in icdict.items():
            # initializes variables as variables in formulas
            self.st.variables[k] = v

    def datasources(self, icdict):
        # key (label) value (source file)
        for k, v in icdict.items():
            data_np_array = np.genfromtxt("user_provided_data/" + v)

            self.data[k] = v


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

