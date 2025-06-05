import numpy as np
from tqdm import tqdm
from PV_Circuit_Model.measurement import *
from matplotlib import pyplot as plt
import copy
import inspect

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

class Fit_Dashboard():
    def __init__(self,nrows,ncols,save_file_name=None,measurements=None,RMS_errors=None):
        _, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6, 5))
        self.axs = axs
        self.nrows = nrows
        self.ncols = ncols
        for ax in axs.flatten():
            ax.set_visible(False)
        self.plot_what = []
        self.define_plot_what(which_axs=0,plot_type="error")
        self.measurements = measurements # pointer
        self.RMS_errors = RMS_errors # pointer
        self.save_file_name = save_file_name
    def define_plot_what(self, which_axs=None, measurement_type=None, key_parameter=None, measurement_condition={}, 
                         plot_type="exp_vs_sim", x_axis=None, title=None, plot_style_parameters={}):
        # plot type is "error", "exp_vs_sim", "overlap_curves" or "overlap_key_parameter"
        # if measurement_type is None, then plot_type defaults to "error"
        # if plot type is "exp_vs_sim" or "overlap_key_parameter", need to specify key_parameter
        # if plot type is "overlap_key_parameter", need to additionally specify x_axis
        if which_axs==None:
            if len(self.plot_what)==0:
                which_axs=0
            else:
                which_axs = self.plot_what[-1]["which_axs"] + 1
        self.plot_what.append({"which_axs":which_axs,
                               "measurement_type":measurement_type, 
                             "key_parameter": key_parameter, 
                             "measurement_condition": measurement_condition,
                             "plot_type": plot_type,
                             "x_axis": x_axis,
                             "title": title,
                             "plot_style_parameters": plot_style_parameters})
    @staticmethod 
    def convert_scatter_valid_kwargs(plot_style_parameters):
        scatter_args = inspect.signature(plt.Axes.scatter).parameters
        kwargs = {}
        for key, value in plot_style_parameters.items():
            if key in scatter_args:
                kwargs[key] = value
        return kwargs
    def plot(self):
        for ax in self.axs.flatten():
            ax.clear()
        for i in range(len(self.plot_what)):
            which_axs = self.plot_what[i]["which_axs"]
            title = self.plot_what[i]["title"]
            ax = self.axs.flatten()[which_axs]
            ax.tick_params(labelsize=6) 
            kwargs = self.convert_scatter_valid_kwargs(self.plot_what[i]["plot_style_parameters"])
            if self.plot_what[i]["plot_type"] == "error":
                ax.set_visible(True)
                ax.scatter(np.arange(0,len(self.RMS_errors)), np.log10(np.array(self.RMS_errors)),s=3,**kwargs) 
                ax.set_title("RMS_error", fontsize=6)
                ax.set_xlabel("Iteration", fontsize=6)
                ax.set_ylabel("log10(Error)", fontsize=6)
            elif self.plot_what[i]["plot_type"] == "overlap_curves":
                measurement_type = self.plot_what[i]["measurement_type"]
                measurement_condition = self.plot_what[i]["measurement_condition"]
                for measurement in self.measurements:
                    if isinstance(measurement,measurement_type):
                        meets_all_conditions = True
                        for key, value in measurement_condition.items():
                            if not (key in measurement.measurement_condition and measurement.measurement_condition[key]==value):
                                meets_all_conditions = False
                        if meets_all_conditions:
                            ax.set_visible(True)
                            measurement_type.plot_func(measurement.measurement_data,color="gray",ax=ax,title="Suns-Voc",kwargs=None)
                            measurement_type.plot_func(measurement.simulated_data,ax=ax,title="Suns-Voc",kwargs=self.plot_what[i]["plot_style_parameters"])
                            if title is not None:
                                ax.set_title(title, fontsize=6)
            else:
                measurement_type = self.plot_what[i]["measurement_type"]
                key_parameter = self.plot_what[i]["key_parameter"]
                measurement_condition = self.plot_what[i]["measurement_condition"]
                categories = []
                values = []
                for key, value in measurement_condition.items():
                    categories.append(key)
                    values.append(value)
                cond_key = self.plot_what[i]["x_axis"]
                result = get_measurements_groups(self.measurements,
                                measurement_class=measurement_type,
                                categories=categories,
                                optional_x_axis=cond_key)
                exp_groups = result[0]
                sim_groups = result[1]
                index = (key_parameter, *values)
                if index in exp_groups:
                    ax.set_visible(True)
                    exp_data = exp_groups[index]
                    sim_data = sim_groups[index]
                    if self.plot_what[i]["x_axis"] is not None:
                        x_axis_groups = result[2]
                        x_axis_data = x_axis_groups[index]
                    match self.plot_what[i]["plot_type"]:
                        case "exp_vs_sim":
                            ax.plot(exp_data,exp_data,color="gray",linewidth=0.5)
                            ax.scatter(exp_data, sim_data,s=3,**kwargs)
                            ax.set_xlabel(key_parameter+"(exp)", fontsize=6)
                            ax.set_ylabel(key_parameter+"(sim)", fontsize=6)
                        case "overlap_key_parameter":
                            ax.scatter(x_axis_data,exp_data,color="gray",s=3)
                            ax.scatter(x_axis_data,sim_data,s=3,**kwargs)
                            ax.set_xlabel(cond_key, fontsize=6)
                            ax.set_ylabel(key_parameter, fontsize=6)
                    if title is not None:
                        ax.set_title(title, fontsize=6)
        plt.tight_layout()
        if plt.gcf().canvas.manager is None:
            plt.show(block=False)
        plt.draw()
        if self.save_file_name is not None:
            word = self.save_file_name + "calibration_fit_round_"+str(len(self.RMS_errors)-1)+".jpg"
            plt.savefig(word, format='jpg', dpi=300)
        plt.pause(0.1)

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
    print("kaka")
    print(fit_parameters.get("value"))
    fit_parameters.set("value",new_values)
    print(fit_parameters.get("value"))
    return new_values

# measurement_samples = collection of devices (Cell, Module, etc)
# each with its measurements stored inside .measurements attribute
# could be one sample only
# could be mulitple samples
def fit_routine(measurement_samples,fit_parameters,
                routine_functions,fit_dashboard=None,
                aux={},num_of_epochs=10):
    # Initial Guess
    routine_functions["initial_guess"](fit_parameters,measurement_samples,aux)
    RMS_errors = []
    record = []
    if fit_dashboard is None:
        fit_dashboard = Fit_Dashboard(2,2)
    fit_dashboard.RMS_errors = RMS_errors
    fit_dashboard.measurements = collate_device_measurements(measurement_samples)
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
                fit_dashboard.plot()
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