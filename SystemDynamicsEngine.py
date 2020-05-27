import networkx as nx
import json
import msgfy
import cexprtk
import numpy as np
import warnings
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
                self.G.add_node(node["key"], stock=True, formula=node["formula"], negative=node["negative"], label=node['label'])
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
                self.file_handlers[node['key']] = np.genfromtxt("user_provided_data/" + node['src']+".csv")
                self.st.variables[node['key']] = 0
                self.G.add_node(node["key"], data=True, source_file=node)
            elif node['category'] == "OfNodes":
                # groups have no semantic meaning when simulating
                pass
            else:
                raise Exception("Undefined node type '" + node['category'] + "' in node " + str(node['label']))
        # create flows
        for flow in js["linkDataArray"]:
            if flow['category'] == 'flow':
                self.G.add_edge(flow["from"], flow["labelKeys"][0])
                self.G.add_edge(flow["labelKeys"][0], flow["to"])
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
        return result_array


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

                # preceeding stock is allowed to go negative, just add the inflow
                if self.G.nodes[predecessor_stock_key]["negative"] == True:
                    inflow_value += predecessors_inflow
                else:
                # preceeding stock is not allowed to go negative
                # check if stock passed 0
                    if predecessor_stock_value[0] == True:
                        inflow_value += 0
                    elif predecessor_stock_value[1][-1] >= predecessors_inflow:
                        inflow_value += predecessors_inflow
                    elif predecessor_stock_value[1][-1] < predecessors_inflow:
                        inflow_value += predecessors_inflow

            # accumulate outflow
            for s in successors:
                outflow_value += cexprtk.Expression(self.G.nodes[s]["formula"], self.st).value()

            # check if stock is allowed to be negative
            net_flow = result_values[i-1] + inflow_value*delta_t - outflow_value*delta_t
            if self.G.nodes[stock_key]["negative"]:
                result_values.append(net_flow)
            else:
                if net_flow > 0:
                    result_values.append(net_flow)
                else:
                    return (True, result_values)

        return (False, result_values)

    # returns a array of simulated values for a certain stock over a certain period
    def eval_stock(self, keys, lower_bound=0, upper_bound=10, delta_t=1, precision=2):
        simulation_results = []
        for k in keys:
            flow = self.compute_flow_analytically(k, t_0=lower_bound, upper_bound=upper_bound, delta_t=delta_t)
            simulation_results.append({"label": k, "values": list(np.round(flow, precision))})
        return simulation_results

    def get_stocks_str(self):
        return list(nx.get_node_attributes(self.G, 'stock').keys())



