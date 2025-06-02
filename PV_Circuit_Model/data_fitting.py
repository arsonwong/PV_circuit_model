import numpy as np
from tqdm import tqdm
from PV_Circuit_Model.measurement import *
from matplotlib import pyplot as plt
import copy

class Fit_Parameter():
    def __init__(self,name="variable",value=0.0,nominal_value=None,d_value=None,abs_min=-np.inf,abs_max=np.inf,is_log=False):
        self.name = name
        self.value = value
        if nominal_value is None:
            self.nominal_value = value
        else:
            self.nominal_value = nominal_value
        self.is_log = is_log
        self.abs_min = abs_min
        self.abs_max = abs_max
        self.set_d_value(d_value)
        self.this_min = -np.inf
        self.this_max = np.inf
        self.enabled = True
        self.is_differential = False
    def initialize(self,value):
        self.nominal_value = value
        self.value = value
    def set_nominal(self):
        self.nominal_value = self.value
    def get_parameter(self):
        value_ = self.value
        if self.is_differential:
            value_ += self.d_value
        if self.is_log:
            return 10**(value_)
        else:
            return value_
    def set_d_value(self,d_value=None):
        if d_value is not None:
            self.d_value = d_value
        else:
            if self.is_log:
                self.d_value = np.log10(2)
            else:
                self.d_value = self.value / 100
    def limit_order_of_mag(self,order_of_mag=1.0):
        if self.is_log:
            self.this_min = self.value - order_of_mag
            self.this_max = self.value + order_of_mag
        else:
            self.this_min = self.value/10**(order_of_mag)
            self.this_max = self.value*10**(order_of_mag)
    def limit_delta(self,delta):
        self.this_min = self.value - delta
        self.this_max = self.value + delta
    def get_min(self):
        return max(self.this_min,self.abs_min)       
    def get_max(self):
        return min(self.this_max,self.abs_max)    
    def check_max_min(self):
        self.nominal_value = max(self.nominal_value,self.abs_min)    
        self.nominal_value = min(self.nominal_value,self.abs_max)  
        self.value = max(self.value,self.abs_min)    
        self.value = min(self.value,self.abs_max)  

class Fit_Parameters():
    def __init__(self,fit_parameters=None,names=None):
        if fit_parameters is not None:
            self.fit_parameters = fit_parameters
        elif names is not None:
            self.fit_parameters = [Fit_Parameter(name=name) for name in names]
        else:
            self.fit_parameters = []
        self.is_differential = False
    def add_fit_parameter(self,fit_parameter):
        self.fit_parameters.append(fit_parameter)
    def enable_parameter(self,name=None):
        for element in self.fit_parameters:
            if name is None or element.name==name:
                element.enabled = False
    def disable_parameter(self,name=None):
        for element in self.fit_parameters:
            if name is None or element.name==name:
                element.enabled = False
    def delete_fit_parameter(self,name):
        for element in self.fit_parameters:
            if element.name==name:
                self.fit_parameters.remove(element)
                break
    def get(self, attribute, names=None, enabled_only=True):
        if names is not None:
            if not isinstance(names,list):
                names = [names]
        list_ = []
        for element in self.fit_parameters:
            if (names is not None and element.name in names) or (names is None and ((not enabled_only) or element.enabled)):
                if attribute=="min":
                    list_.append(element.get_min())
                elif attribute=="max":
                    list_.append(element.get_max())
                else:
                    list_.append(getattr(element, attribute))
        if len(list_)==1:
            return list_[0]
        return list_
    def set(self, attribute, values, names=None, enabled_only=True):
        if names is not None:
            if not isinstance(names,list):
                names = [names]
        if not isinstance(values,list) and not isinstance(values,np.ndarray):
            values = [values]*self.num_of_parameters()
        count = 0
        for element in self.fit_parameters:
            if (names is not None and element.name in names) or (names is None and ((not enabled_only) or element.enabled)):
                setattr(element, attribute,values[count])
                element.check_max_min()
                count += 1
    def initialize(self, values, names=None, enabled_only=True):
        self.set("value",values,names=names,enabled_only=enabled_only)
        self.set("nominal_value",values,names=names,enabled_only=enabled_only)
    def set_nominal(self):
        for element in self.fit_parameters:
            element.set_nominal()
    def set_d_value(self):
        for element in self.fit_parameters:
            element.set_d_value()
    def set_differential(self,which=-1,enabled_only=True):
        count = 0
        if which<0:
            self.is_differential = False
        else:
            self.is_differential = True
        for element in self.fit_parameters:
            element.is_differential = False
            if (not enabled_only) or element.enabled:
                if count==which:
                    element.is_differential = True
                count += 1
    def get_parameters(self):
        dict_ = {}
        for element in self.fit_parameters:
            dict_[element.name] = element.get_parameter()
        return dict_
    def limit_order_of_mag(self,order_of_mag=1.0):
        for element in self.fit_parameters:
            element.limit_order_of_mag(order_of_mag=order_of_mag)
    def num_of_parameters(self):
        return len(self.fit_parameters)
    def num_of_enabled_parameters(self):
        count = 0
        for element in self.fit_parameters:
            if element.enabled:
                count += 1
        return count
    def apply_to_device(self, device):
        pass
    def __str__(self):
        return str(self.get_parameters())

def linear_regression(M, Y, fit_parameters, aux={}): 
    alpha = 1e-5
    regularization_method=0 
    if "alpha" in aux:
        alpha = aux["alpha"]
    if "regularization_method" in aux:
        regularization_method = aux["regularization_method"]
    if "limit_order_of_mag" in aux:
        if aux["limit_order_of_mag"]:
            fit_parameters.limit_order_of_mag()
    min_values = fit_parameters.get("min")
    max_values = fit_parameters.get("max")
    values = np.array(fit_parameters.get("value"))
    nominal_values = np.array(fit_parameters.get("nominal_value"))
    dvalues = np.array(fit_parameters.get("d_value"))
    too_high_indices = []
    too_low_indices = []  
    included = np.ones_like(nominal_values)
    included_indices = np.where(included==1)[0]
    Ybias = np.zeros_like(Y)
    Xbias = np.zeros_like(nominal_values)

    while True:
        Y_ = Y - Ybias
        M_ = M[:,included_indices]
        Y2 = np.vstack([Y_[:,None],
                        (alpha*(nominal_values[included_indices] - values[included_indices])/dvalues[included_indices])[:,None]])
        M2 = np.vstack([M_,alpha*np.identity(M_.shape[1])])
        # another regularization on how much variables can change at a time
        alpha2 = 1e-7
        memory = [[too_high_indices.copy(),too_low_indices.copy()],included.copy(),included_indices.copy(),Ybias.copy(),Xbias.copy(),M2.copy(),Y2.copy()]
        len_excluded_indices = len(too_low_indices)+len(too_high_indices)
        while True:
            too_high_indices = memory[0][0].copy()
            too_low_indices = memory[0][1].copy()
            included = memory[1].copy()
            included_indices = memory[2].copy()
            Ybias = memory[3].copy()   
            Xbias = memory[4].copy()
            M2 = memory[5].copy()
            Y2 = memory[6].copy()
            M2 = np.vstack([M_,alpha2*np.identity(M_.shape[1])])
            Y2 = np.vstack([Y_[:,None],np.zeros((M_.shape[1],1))])

            MTM = M2.T @ M2
            MTY = M2.T @ Y2
            X_ = np.linalg.solve(MTM, MTY)

            X = Xbias.copy()
            X[included_indices] = X_[:,0]

            delta = X*dvalues
            new_values = values + delta

            find_ = np.where(new_values < min_values)[0]
            if len(find_) > 0:
                too_low_indices.extend(find_)
            find2_ = np.where(new_values > max_values)[0]
            if len(find2_) > 0:
                too_high_indices.extend(find2_)
            if len(too_low_indices) > 0:
                too_low_indices = list(np.unique(np.array(too_low_indices)))

            if regularization_method==1 or len(too_low_indices)+len(too_high_indices) == len_excluded_indices:
                break
            alpha2 *= 3
        
        if len(too_low_indices)+len(too_high_indices) == len_excluded_indices:
            break
        Xbias[too_low_indices] = (min_values[too_low_indices]-values[too_low_indices])/dvalues[too_low_indices]
        Xbias[too_high_indices] = (max_values[too_high_indices]-values[too_high_indices])/dvalues[too_high_indices]
        Ybias = M @ Xbias
        included = np.ones_like(values)
        included[too_low_indices] = 0
        included[too_high_indices] = 0
        included_indices = np.where(included==1)[0]
        assert(len(included_indices)>0)
    fit_parameters.set("value",new_values)
    return new_values

def plot_error(aux,output,RMS_errors):
    if "axs" not in aux:
        measurements = collate_device_measurements(sample)
        _, aux["axs"] = plt.subplots(nrows=2, ncols=2, figsize=(6, 5))
        for ax in aux["axs"].flatten():
            ax.set_visible(False)
    if "axs" in aux:
        axs = aux["axs"]
        ax = axs.flatten()[0]
        ax.set_visible(True)
        ax.clear()
        ax.scatter(np.arange(0,len(RMS_errors)), np.log10(np.array(RMS_errors)),s=3)
        ax.tick_params(labelsize=6)  
        ax.set_title("Error", fontsize=6)
        ax.set_xlabel("Iteration", fontsize=6)
        ax.set_ylabel("log10(Error)", fontsize=6)
        ax = axs.flatten()[1]
        ax.set_visible(True)
        ax.clear()
        
        plt.tight_layout()
        if plt.gcf().canvas.manager is None:
            plt.show(block=False)
        plt.draw()
        if "prefix" in aux:
            word = aux["prefix"] + "calibration_fit_round_"+str(len(RMS_errors)-1)+".jpg"
            plt.savefig(word, format='jpg', dpi=300)
        plt.pause(0.1)

# measurement_samples = collection of devices (Cell, Module, etc)
# each with its measurements stored inside .measurements attribute
# could be one sample only
# could be mulitple samples
def fit_routine(measurement_samples,fit_parameters,
                routine_functions,
                aux={},num_of_epochs=10):
    # Initial Guess
    routine_functions["initial_guess"](fit_parameters,measurement_samples,aux)
    RMS_errors = []
    record = []
    if "comparison_function_iterations" not in aux:
        aux["comparison_function_iterations"] = 1
    aux["pbar"] = tqdm(total=((num_of_epochs-1)*(fit_parameters.num_of_enabled_parameters()+1)+1)*aux["comparison_function_iterations"],desc="Calibrating")
    for epoch in range(num_of_epochs):
        M = []
        for iteration in range(fit_parameters.num_of_enabled_parameters()+1):
            fit_parameters.set_differential(iteration-1)
            pbar_before = aux["pbar"].n
            output = routine_functions["comparison_function"](fit_parameters,measurement_samples,aux)
            pbar_after = aux["pbar"].n
            if pbar_after == pbar_before:
                aux["pbar"].update(aux["comparison_function_iterations"])
            if iteration==0:
                Y = np.array(output["error_vector"])
                RMS_errors.append(np.sqrt(np.mean(Y**2)))
                record.append({"fit_parameters": copy.deepcopy(fit_parameters),"output": output})
                if "plot_function" in routine_functions:
                    routine_functions["plot_function"](aux,output,RMS_errors)
            else:
                M.append(output["differential_vector"])
            if epoch==num_of_epochs-1:
                index = np.argmin(np.array(RMS_errors))
                fit_parameters = record[index]["fit_parameters"]
                output = record[index]["output"]
                aux["pbar"].close()
                return output
        M = np.array(M)
        M = M.T
        fit_parameters.set_differential(-1)
        routine_functions["update_function"](M, Y, fit_parameters, aux)
        if "post_update_function" in routine_functions:
            routine_functions["post_update_function"](fit_parameters)