import numpy as np
from matplotlib import pyplot as plt
from PV_Circuit_Model.cell_analysis import *
from PV_Circuit_Model.cell import *
from PV_Circuit_Model.multi_junction_cell import *
import numbers

class Measurement():
    keys = []
    # measurement can be on its own, or belonging to a device
    def __init__(self,measurement_condition,measurement_data,device=None):
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

def get_measurements_groups(measurements,measurement_class=None,
                            include_tags=None,exclude_tags=None,categories=[],
                            optional_x_axis=None):
    exp_groups = {}
    sim_groups = {}
    x_axis_groups = {}
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
                    sim_ = []
                    if key in measurement.simulated_key_parameters:
                        sim_ = measurement.simulated_key_parameters[key]
                    if optional_x_axis is not None:
                        condition = measurement.measurement_condition[optional_x_axis]
                    if isinstance(exp_,numbers.Number):
                        exp_ = [exp_]
                    elif isinstance(exp_, np.ndarray):
                        exp_ = exp_.tolist()
                    if isinstance(sim_,numbers.Number):
                        sim_ = [sim_]
                    elif isinstance(sim_, np.ndarray):
                        sim_ = sim_.tolist()
                    exp_groups[tuple_].extend(exp_)
                    sim_groups[tuple_].extend(sim_)
                    if optional_x_axis is not None:
                        x_axis_groups[tuple_].extend([condition]*len(exp_))
    if optional_x_axis is not None:
        return exp_groups, sim_groups, x_axis_groups
    return exp_groups,sim_groups

def collate_device_measurements(devices,measurement_class=None,include_tags=None,exclude_tags=None):
    measurement_list = []
    if not isinstance(devices,list):
        devices = [devices]
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
    keys = ["Voc", "Isc", "Pmax"]
    def __init__(self,Suns,IV_curve,is_dark=False,temperature=25,IL=None,JL=None,**kwargs):
        super().__init__(measurement_condition={'Suns':Suns,'IL':IL,'JL':JL,'is_dark':is_dark,
                                                'temperature':temperature},
                         measurement_data=IV_curve,**kwargs)
    @staticmethod
    def derive_key_parameters(data,key_parameters,conditions):
        if not conditions["is_dark"]:
            key_parameters["Voc"] = get_Voc(data)
            key_parameters["Isc"] = get_Isc(data)
            key_parameters["Pmax"], _, _ = get_Pmax(data)
        else:
            Rshunt = Rshunt_extraction(data)
            key_parameters["log_shunt_cond"] = np.log10(1/Rshunt)
    def simulate(self,device=None):
        temperature = self.measurement_condition["temperature"]
        Suns = self.measurement_condition["Suns"]
        IL = self.measurement_condition["IL"]
        JL = self.measurement_condition["JL"]
        if device is None:
            device = self.parent_device
        device.set_temperature(temperature,rebuild_IV=False)
        if JL is not None:
            device.set_JL(JL,rebuild_IV=False)
            device.set_Suns(1.0)
        elif IL is not None:
            device.set_IL(IL,rebuild_IV=False)
            device.set_Suns(1.0)
        else:
            device.set_Suns(Suns)
        self.simulated_data = device.IV_table
        self.derive_key_parameters(self.simulated_data, self.simulated_key_parameters, self.measurement_condition)
    @staticmethod
    def plot_func(data,color="black"):
        plt.plot(data[0,:],data[1,:],color=color)
        _, Vmp, Imp = get_Pmax(data)
        plt.scatter(Vmp,Imp)
        plt.xlabel("Voltage (V)")
        plt.ylabel("Current (A)")

class dark_IV_measurement(IV_measurement):
    keys = ["log_shunt_cond"]
    @staticmethod
    def derive_key_parameters(data,key_parameters,conditions):
        Rshunt = Rshunt_extraction(data)
        key_parameters["log_shunt_cond"] = np.log10(1/Rshunt)

class Suns_Voc_measurement(Measurement):
    keys = ["Voc"]
    def __init__(self,Suns_Isc_Voc_curve,temperature=25,**kwargs):
        super().__init__(measurement_condition={'temperature':temperature},
                         measurement_data=Suns_Isc_Voc_curve,**kwargs)
    @staticmethod
    def derive_key_parameters(data,key_parameters,conditions):
        key_parameters["Voc"] = data[:,0]
    def simulate(self,device=None):
        num_col = self.measurement_data.shape[1]
        num_subcells = int((num_col-1)/2)
        Suns = self.measurement_data[:,1:num_subcells+1]
        Iscs = self.measurement_data[:,num_subcells+1:]
        if np.isnan(Suns[0,0]):
            Suns = None
        if np.isnan(Iscs[0,0]):
            Iscs = None
        if device is None:
            device = self.parent_device
        self.simulated_data, _ = simulate_Suns_Voc(device, Suns=Suns, Iscs=Iscs)
        self.derive_key_parameters(self.simulated_data, self.simulated_key_parameters, None)
    @staticmethod
    def plot_func(data,color="black"):
        num_col = data.shape[1]
        num_subcells = int((num_col-1)/2)
        y_label = "log10(Suns)"
        ys = np.max(data[:,1:num_subcells+1],axis=1)
        if np.isnan(ys[0]):
            y_label = "log10(Current(A))"
            ys = np.max(data[:,num_subcells+1:],axis=1)
        ys = np.log10(ys)
        plt.scatter(data[:,0],ys,color=color)
        plt.xlabel("Voc (V)")
        plt.ylabel(y_label)

def simulate_Suns_Voc(cell, Suns=None, Iscs=None):
    subcells_num = 1
    if isinstance(cell,MultiJunctionCell):
        subcells_num = len(cell.cells)
    if Suns is None and Iscs is None:
        Suns = 10.0**(np.arange(-3,1,0.1))
        Suns = Suns[:,None]
    if Iscs is not None:
        if isinstance(Iscs,numbers.Number):
            Iscs = Iscs*np.ones((1,subcells_num))
        if Iscs.ndim == 1:
            Iscs = Iscs[:,None]
        assert(Iscs.shape[1]==subcells_num)
        Suns = np.ones_like(Iscs)*np.NaN
    else:
        if isinstance(Suns,numbers.Number):
            Suns = Suns*np.ones((1,subcells_num))
        if Suns.ndim == 1:
            Suns = Suns[:,None]
        assert(Suns.shape[1]==subcells_num)
        Iscs = np.ones_like(Suns)*np.NaN
    Vocs = []
    cell.set_Suns(1.0, rebuild_IV=False)
    for i, _ in enumerate(Suns[:,0]):
        if not np.isnan(Suns[i,0]):
            if isinstance(cell,MultiJunctionCell):
                for j, cell_ in enumerate(cell.cells):
                    cell_.set_Suns(Suns[i,j])
                cell.build_IV()
            else:
                cell.set_Suns(Suns[i,0])
        else:
            if isinstance(cell,MultiJunctionCell):
                for j, cell_ in enumerate(cell.cells):
                    cell_.set_IL(Iscs[i,j], temperature=cell.temperature)
                cell.build_IV()
                # for cell_ in cell.cells:
                #     cell_.plot()
                # cell.plot()
                # cell.show()
            else:
                cell.set_IL(Iscs[i,0], temperature=cell.temperature)
        Vocs.append(cell.get_Voc())
    Suns_Isc_Voc_curve = np.hstack([np.array(Vocs)[:,None],Suns,Iscs])
    if len(Vocs)==1:
        Vocs = Vocs[0]
    return Suns_Isc_Voc_curve, Vocs


