import cexprtk
import numpy as np
from scipy.integrate import odeint
import numpy.polynomial.polynomial as polynomial

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

        def ite(a, b, c):
            return b if a else c
        self.st.functions['ite'] = ite


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





