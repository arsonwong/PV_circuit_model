import numpy as np
from PV_Circuit_Model.cell_analysis import *
import numbers

class Measurement():
    # measurement can be on its own, or belonging to a device
    def __init__(self,measurement_condition,measurement_data,device=None):
        if not hasattr(self,"keys"):
            self.keys = []
        self.measurement_condition = measurement_condition
        self.measurement_data = measurement_data
        self.simulated_data = None
        self.tag = None
        self.key_parameters = {}
        self.simulated_key_parameters = {}
        self.simulated_key_parameters_baseline = {}
        self.unit_errors = {}
        self.parent_device=device
        self.derive_key_parameters(self.measurement_data,self.key_parameters,self.measurement_condition)
        self.set_unit_errors()
    def set_unit_error(self,key,value=None):
        if value is not None:
            self.unit_errors[key] = value
        else: # just make as percentage
            if isinstance(self.key_parameters[key],numbers.Number):
                avg_ = abs(self.key_parameters[key])
                self.unit_errors[key] = avg_
            else: 
                avg_ = np.mean(np.abs(self.key_parameters[key])) 
                self.unit_errors[key] = avg_*np.sqrt(self.key_parameters[key].size)
    def set_unit_errors(self):
        for key in self.keys:
            self.set_unit_error(key)
    @staticmethod
    def derive_key_parameters(data, key_parameters, conditions):
        pass
    def simulate(self,device=None):
        pass
    def plot(self):
        self.plot_func(self.measurement_data,color="blue")
        if self.simulated_data is not None:
            self.plot_func(self.simulated_data,color="red")
        plt.show()
    @staticmethod
    def plot_func(data):
        pass
    def get_diff_vector(self,key_parameters1,key_parameters2):
        diff = []
        for key in self.keys:
            parameter1 = key_parameters1[key]
            parameter2 = key_parameters2[key]
            if isinstance(parameter1,numbers.Number):
                parameter1 = [parameter1]
                parameter2 = [parameter2]
            if isinstance(parameter1, list):
                parameter1 = np.array(parameter1)
                parameter2 = np.array(parameter2)
            has_nan = np.isnan(parameter1).any()
            if has_nan:
                print(type(self))
                print(self.key_parameters)
                print(self.simulated_data)
                assert(1==0)
            has_nan = np.isnan(parameter2).any()
            if has_nan:
                print(type(self))
                print(self.key_parameters)
                print(self.simulated_data)
                assert(1==0)
            diff_vector = parameter1-parameter2
            diff_vector = diff_vector.tolist()
            diff.extend(diff_vector/self.unit_errors[key])
        return np.array(diff)
    def set_simulation_baseline(self):
        self.simulated_key_parameters_baseline = self.simulated_key_parameters.copy()
    def get_error_vector(self):
        return self.get_diff_vector(self.key_parameters,self.simulated_key_parameters)
    def get_differential_vector(self):
        return self.get_diff_vector(self.simulated_key_parameters,self.simulated_key_parameters_baseline)
    def __str__(self):
        return str(self.key_parameters)
    
def assign_measurements(self:CircuitGroup, measurements):
    for measurement in measurements:
        measurement.parent_device = self
    self.measurements = measurements
CircuitGroup.assign_measurements = assign_measurements

def set_simulation_baseline(measurements):
    for measurement in measurements:
        measurement.set_simulation_baseline()
        
def get_measurements_error_vector(measurements,measurement_class=None,include_tags=None,exclude_tags=None):
    vector = []
    for measurement in measurements:
        if measurement_class is None or isinstance(measurement,measurement_class):
            if (include_tags==None or (measurement.tag is not None and measurement.tag in include_tags)) and (exclude_tags==None or (measurement.tag is None or measurement.tag not in exclude_tags)):
                vector.extend(measurement.get_error_vector())
    return vector

def get_measurements_differential_vector(measurements,measurement_class=None,include_tags=None,exclude_tags=None):
    vector = []
    for measurement in measurements:
        if measurement_class is None or isinstance(measurement,measurement_class):
            if (include_tags==None or (measurement.tag is not None and measurement.tag in include_tags)) and (exclude_tags==None or (measurement.tag is None or measurement.tag not in exclude_tags)):
                vector.extend(measurement.get_differential_vector())
    return vector

def get_measurements_groups(measurements,measurement_class=None,include_tags=None,exclude_tags=None,categories=[]):
    exp_groups = {}
    sim_groups = {}
    for measurement in measurements:
        if measurement_class is None or isinstance(measurement,measurement_class):
            if (include_tags==None or (measurement.tag is not None and measurement.tag in include_tags)) and (exclude_tags==None or (measurement.tag is None or measurement.tag not in exclude_tags)):
                conditions = []
                for sought_category in categories:
                    for category, condition in measurement.measurement_condition.items():
                        if category==sought_category:
                            conditions.append(condition)
                for key in measurement.keys:
                    tuple_ = key
                    if len(conditions)>0:
                        tuple_ = tuple([key]+conditions)
                    if tuple_ not in exp_groups:
                        exp_groups[tuple_] = []
                        sim_groups[tuple_] = []
                    exp_ = measurement.key_parameters[key]
                    sim_ = measurement.simulated_key_parameters[key]
                    if isinstance(exp_,numbers.Number):
                        exp_ = [exp_]
                        sim_ = [sim_]
                    elif isinstance(exp_, np.ndarray):
                        exp_ = exp_.tolist()
                        sim_ = sim_.tolist()
                    exp_groups[tuple_].extend(exp_)
                    sim_groups[tuple_].extend(sim_)
    return exp_groups,sim_groups

def collate_device_measurements(devices,measurement_class=None,include_tags=None,exclude_tags=None):
    measurement_list = []
    for device in devices:
        measurements = device.measurements
        for measurement in measurements:
            if measurement_class is None or isinstance(measurement,measurement_class):
                if (include_tags==None or (measurement.tag is not None and measurement.tag in include_tags)) and (exclude_tags==None or (measurement.tag is None or measurement.tag not in exclude_tags)):
                    measurement_list.append(measurement)
    return measurement_list

def simulate_device_measurements(devices,measurement_class=None,include_tags=None,exclude_tags=None):
    for device in devices:
        measurements = device.measurements
        for measurement in measurements:
            if measurement_class is None or isinstance(measurement,measurement_class):
                if (include_tags==None or (measurement.tag is not None and measurement.tag in include_tags)) and (exclude_tags==None or (measurement.tag is None or measurement.tag not in exclude_tags)):
                    measurement.simulate(device)

class IV_measurement(Measurement):
    def __init__(self,Suns,IV_curve,temperature=25,**kwargs):
        self.keys = ["Voc", "Isc", "Pmax"]
        super().__init__(measurement_condition={'Suns':Suns,'temperature':temperature},
                         measurement_data=IV_curve,**kwargs)
    @staticmethod
    def derive_key_parameters(data,key_parameters,conditions):
        key_parameters["Voc"] = get_Voc(data)
        key_parameters["Isc"] = get_Isc(data)
        key_parameters["Pmax"], _, _ = get_Pmax(data)
    def simulate(self,device=None):
        temperature = self.measurement_condition["temperature"]
        Suns = self.measurement_condition["Suns"]
        if device is None:
            device = self.parent_device
        device.set_temperature(temperature,rebuild_IV=False)
        device.set_Suns(Suns)
        self.simulated_data = device.IV_table
        self.derive_key_parameters(self.simulated_data, self.simulated_key_parameters, None)
    @staticmethod
    def plot_func(data,color="black"):
        plt.plot(data[0,:],data[1,:],color=color)
        _, Vmp, Imp = get_Pmax(data)
        plt.scatter(Vmp,Imp)
        plt.xlabel("Voltage (V)")
        plt.ylabel("Current (A)")

class Suns_Voc_measurement(Measurement):
    def __init__(self,Suns_Isc_Voc_curve,temperature=25,**kwargs):
        self.keys = ["Voc"]
        super().__init__(measurement_condition={'temperature':temperature},
                         measurement_data=Suns_Isc_Voc_curve,**kwargs)
    @staticmethod
    def derive_key_parameters(data,key_parameters,conditions):
        key_parameters["Voc"] = data[-1,:]
    def simulate(self,device=None):
        Suns = self.measurement_data[0,:]
        Iscs = self.measurement_data[1,:]
        if np.isnan(Suns[0]):
            Suns = None
        if np.isnan(Iscs[0]):
            Iscs = None
        if device is None:
            device = self.parent_device
        self.simulated_data = simulate_Suns_Voc(device, Suns=Suns, Iscs=Iscs)
        self.derive_key_parameters(self.simulated_data, self.simulated_key_parameters, None)
    @staticmethod
    def plot_func(data):
        y_label = "log10(Suns)"
        ys = data[0,:]
        if np.isnan(ys[0]):
            y_label = "log10(Current(A))"
            ys = data[1,:]
        ys = np.log10(ys)
        plt.plot(data[-1,:],ys)
        plt.xlabel("Voc (V)")
        plt.ylabel(y_label)

def simulate_Suns_Voc(cell:Cell, Suns=None, Iscs=None):
    if Suns is None and Iscs is None:
        Suns = 10.0**(np.arange(-3,1,0.1))
    if Iscs is not None:
        Suns = np.ones_like(Iscs)*np.NaN
    else:
        Iscs = np.ones_like(Suns)*np.NaN
    Vocs = []
    parent_Vocs = [] # if this is a subcell to a tandem cell, report tandem cell as well
    is_in_tandem = False
    if cell.parent is not None and cell.parent.connection=="series":
        is_in_tandem = True
    for i, _ in enumerate(Suns):
        if not np.isnan(Suns[i]):
            cell.set_Suns(Suns[i])
        else:
            cell.set_Suns(1.0)
            cell.set_IL(Iscs[i], temperature=cell.temperature)
        Vocs.append(cell.get_Voc())
        if is_in_tandem:
            cell.parent.build_IV()
            parent_Vocs.append(cell.parent.get_Voc())
    if is_in_tandem:
        Suns_Isc_Voc_curve = np.array([Suns,Iscs,Vocs,parent_Vocs])
    else:
        Suns_Isc_Voc_curve = np.array([Suns,Iscs,Vocs])
    return Suns_Isc_Voc_curve
