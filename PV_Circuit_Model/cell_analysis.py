import numpy as np
from PV_Circuit_Model.circuit_model import *
from PV_Circuit_Model.cell import *
from PV_Circuit_Model.module import *

def get_Voc(argument):
    if isinstance(argument,CircuitGroup):
        IV_curve = argument.IV_table
    else:
        IV_curve = argument
    return np.interp(0,IV_curve[1,:],IV_curve[0,:])
CircuitGroup.get_Voc = get_Voc

def get_Isc(argument):
    if isinstance(argument,CircuitGroup):
        IV_curve = argument.IV_table
    else:
        IV_curve = argument
    return -np.interp(0,IV_curve[0,:],IV_curve[1,:])
CircuitGroup.get_Isc = get_Isc

def get_Pmax(argument, return_op_point=False):
    if isinstance(argument,CircuitGroup):
        IV_curve = argument.IV_table
    else:
        IV_curve = argument
        return_op_point = True
    Voc = get_Voc(IV_curve)
    V = np.linspace(0,Voc,1000)
    I = np.interp(V,IV_curve[0,:],IV_curve[1,:])
    power = -V*I
    index = np.argmax(power)
    Vmp = V[index]
    Imp = I[index]
    max_power = power[index]
    if return_op_point:
        return max_power, Vmp, Imp
    return max_power
CircuitGroup.get_Pmax = get_Pmax

def get_FF(argument):
    if isinstance(argument,CircuitGroup):
        IV_curve = argument.IV_table
    else:
        IV_curve = argument
    Voc = get_Voc(IV_curve)
    Isc = get_Isc(IV_curve)
    Pmax, _, _ = get_Pmax(IV_curve)
    return Pmax/(Isc*Voc)
CircuitGroup.get_FF = get_FF

def Rs_extraction_two_light_IVs(IV_curves):
    Isc0 = -1*get_Isc(IV_curves[0])
    Isc1 = -1*get_Isc(IV_curves[1])
    _, Vmp0, Imp0 = get_Pmax(IV_curves[0])
    delta_I = -Isc0+Imp0
    delta_Is_halfSun = -Isc1+IV_curves[1][1,:]
    V_point = np.interp(delta_I,delta_Is_halfSun,IV_curves[1][0,:])
    Rs = (Vmp0-V_point)/(Isc0-Isc1)
    return Rs

def Rshunt_extraction(IV_curve):
    indices = np.where((IV_curve[0,:]>0.0) & (IV_curve[0,:]<=0.1))[0]
    m, _ = np.polyfit(IV_curve[0,indices], IV_curve[1,indices], deg=1)
    if m <= 0:
        Rshunt = 100000
    else:
        Rshunt = 1/m
    Rshunt = min(Rshunt,100000)
    return Rshunt

def estimate_cell_J01_J02(Jsc,Voc,Pmax=None,FF=1.0,Rs=0.0,Rshunt=1e6,
                          temperature=25,Sun=1.0,thickness=180e-4,Si_intrinsic_limit=True):
    if Pmax is None:
        Pmax = Jsc*Voc*FF          
    VT = get_VT(temperature)
    max_J01 = Jsc/np.exp(Voc/VT)
    for inner_k in range(100):
        trial_cell = make_solar_cell(Jsc, max_J01, 0.0, Rshunt, 
                                     Rs, thickness=thickness,Si_intrinsic_limit=Si_intrinsic_limit)
        trial_cell.set_temperature(temperature,rebuild_IV=False)
        trial_cell.set_Suns(Sun)
        Voc_ = trial_cell.get_Voc()
        if abs(Voc_-Voc) < 1e-4:
            break 
        max_J01 *= np.exp((Voc_-Voc)/VT)
    max_J02 = Jsc/np.exp(Voc/(2*VT))
    for inner_k in range(100):
        trial_cell = make_solar_cell(Jsc, 0.0, max_J02, Rshunt, Rs, 
                                     thickness=thickness,Si_intrinsic_limit=Si_intrinsic_limit)
        trial_cell.set_temperature(temperature,rebuild_IV=False)
        trial_cell.set_Suns(Sun)
        Voc_ = trial_cell.get_Voc()
        if abs(Voc_-Voc) < 1e-4:
            break 
        max_J02 *= np.exp((Voc_-Voc)/(2*VT))
    outer_record = []
    for outer_k in range(100):
        if outer_k==0:
            trial_J01 = 0.0
        elif outer_k==1:
            trial_J01 = max_J01
        else:
            outer_record_ = np.array(outer_record)
            indices = np.argsort(outer_record_[:,0])
            outer_record_ = outer_record_[indices,:]
            trial_J01 = interp_(Pmax, outer_record_[:,1], outer_record_[:,0])[0]
            trial_J01 = max(trial_J01, 0.0)
            trial_J01 = min(trial_J01, max_J01)
        inner_record = []
        for inner_k in range(100):
            if inner_k==0:
                trial_J02 = 0.0
            elif inner_k==1:
                trial_J02 = max_J02
            else:
                inner_record_ = np.array(inner_record)
                indices = np.argsort(inner_record_[:,0])
                inner_record_ = inner_record_[indices,:]
                trial_J02 = interp_(Voc, inner_record_[:,1], inner_record_[:,0])[0]
                trial_J02 = max(trial_J02, 0.0)
                trial_J02 = min(trial_J02, max_J02)
            trial_cell = make_solar_cell(Jsc, trial_J01, trial_J02, Rshunt, Rs,thickness=thickness,
                                         Si_intrinsic_limit=Si_intrinsic_limit)
            trial_cell.set_temperature(temperature,rebuild_IV=False)
            trial_cell.set_Suns(Sun)
            Voc_ = trial_cell.get_Voc()
            if abs(Voc_-Voc) < 1e-4 or (trial_J02==0 and Voc_<Voc) or (trial_J02==max_J02 and Voc_>Voc):
                break 
            inner_record.append([trial_J02,Voc_])
        Pmax_ = trial_cell.get_Pmax()
        outer_record.append([trial_J01,Pmax_])
        if abs(Voc_-Voc)<1e-4 and abs(Pmax_-Pmax)/Pmax<1e-4:
            break
        if outer_k==1 and Pmax_ < Pmax: # will never be bigger then
            break
    return trial_J01, trial_J02

def plot(self, fourth_quadrant=True, show_IV_parameters=True):
    if fourth_quadrant and isinstance(self,CircuitGroup):
        Voc = self.get_Voc()
        Isc = self.get_Isc()
        plt.plot(self.IV_table[0,:],-self.IV_table[1,:])
        if self.operating_point is not None:
            plt.plot(self.operating_point[0],-self.operating_point[1],marker='o')
            if len(self.operating_point)==3:
                plt.plot(self.operating_point[2],-self.operating_point[1],marker='o')
        plt.xlim((0,Voc*1.1))
        plt.ylim((0,Isc*1.1))
    else:
        plt.plot(self.IV_table[0,:],self.IV_table[1,:])
        if self.operating_point is not None:
            plt.plot(self.operating_point[0],self.operating_point[1],marker='o')
            if len(self.operating_point)==3:
                plt.plot(self.operating_point[2],self.operating_point[1],marker='o')
    if show_IV_parameters and fourth_quadrant and (isinstance(self,Cell) or isinstance(self,Module)):
        max_power, Vmp, Imp = self.get_Pmax(return_op_point=True)
        Voc = self.get_Voc()
        Isc = self.get_Isc()
        FF = self.get_FF()
        y_space = 0.07
        plt.plot(Voc,0,marker='o',color="blue")
        plt.plot(0,Isc,marker='o',color="blue")
        if fourth_quadrant:
            Imp *= -1
        plt.plot(Vmp,Imp,marker='o',color="blue")
        if isinstance(self,Cell):
            plt.text(Voc*0.05, Isc*(0.8-0*y_space), f"Isc = {Isc:.3f} A")
            plt.text(Voc*0.05, Isc*(0.8-1*y_space), f"Jsc = {Isc/self.area*1000:.3f} mA/cm2")
            plt.text(Voc*0.05, Isc*(0.8-2*y_space), f"Voc = {Voc:.4f} V")
            plt.text(Voc*0.05, Isc*(0.8-3*y_space), f"FF = {FF*100:.3f} %")
            plt.text(Voc*0.05, Isc*(0.8-4*y_space), f"Pmax = {max_power:.3f} W")
            plt.text(Voc*0.05, Isc*(0.8-5*y_space), f"Eff = {max_power/self.area*1000:.3f} %")
            plt.text(Voc*0.05, Isc*(0.8-6*y_space), f"Area = {self.area:.3f} cm2")
        else:
            plt.text(Voc*0.05, Isc*(0.8-0*y_space), f"Isc = {Isc:.3f} A")
            plt.text(Voc*0.05, Isc*(0.8-1*y_space), f"Voc = {Voc:.2f} V")
            plt.text(Voc*0.05, Isc*(0.8-2*y_space), f"FF = {FF*100:.3f} %")
            plt.text(Voc*0.05, Isc*(0.8-3*y_space), f"Pmax = {max_power:.2f} W")
    plt.xlabel("Voltage (V)")
    plt.ylabel("Current (A)")
CircuitGroup.plot = plot
CircuitElement.plot = plot

def show(self):
    plt.show()
CircuitGroup.show = show