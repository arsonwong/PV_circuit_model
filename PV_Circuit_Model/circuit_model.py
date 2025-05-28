import numpy as np
from matplotlib import pyplot as plt
from PV_Circuit_Model.utilities import *
import copy
from tqdm import tqdm

pbar = None
x_spacing = 1.5
y_spacing = 0.2

class const():
    VT = 0.02568 # thermal voltage at 25C

def get_VT(temperature):
    return const.VT*(temperature + 273.15)/(25 + 273.15)

class CircuitElement:
    def __init__(self,tag=None):
        self.IV_table = None
        self.tag = tag
        self.operating_point = None #V,I
        self.circuit_diagram_extent = [0, 0.8]
        self.parent = None
    def set_operating_point(self,V=None,I=None,dark_IV=None):
        if V is not None:
            I = np.interp(V,self.IV_table[0,:],self.IV_table[1,:])
        elif I is not None:
            V = np.interp(I,self.IV_table[1,:],self.IV_table[0,:])
        self.operating_point = [V,I]
    def get_value_text(self):
        pass
    def get_draw_func(self):
        pass
    def draw(self, ax=None, x=0, y=0, color="black", display_value=False):
        text = None
        if display_value:
            text = self.get_value_text()
        draw_symbol(self.get_draw_func(),ax=ax,x=x,y=y,color=color,text=text)
    def null_IV(self, keep_dark=False):
        self.IV_table = None
        if self.parent is not None:
            self.parent.null_IV(keep_dark=keep_dark)

class CurrentSource(CircuitElement):
    def __init__(self, IL, Suns=1.0, temperature=25, temp_coeff=0.0, tag=None):
        super().__init__(tag=tag)
        if np.isnan(IL):
            assert(1==0)
        self.IL = IL
        self.refSuns = Suns
        self.Suns = Suns
        self.refIL = IL
        self.refT = temperature
        self.T = temperature
        self.temp_coeff = temp_coeff

    def set_IL(self,IL):
        self.IL = IL
        self.null_IV(keep_dark=True)

    def copy(self,source):
        self.refSuns = source.refSuns
        self.Suns = source.Suns
        self.refIL = source.refIL
        self.refT = source.refT
        self.T = source.T
        self.temp_coeff = source.temp_coeff
        self.set_IL(source.IL)

    def changeTemperatureAndSuns(self,temperature=None,Suns=None,rebuild_IV=True):
        if Suns is not None:
            self.Suns = Suns
        if temperature is not None:
            self.T = temperature
        self.set_IL(self.Suns*(self.refIL / self.refSuns + self.temp_coeff * (self.T - self.refT)))
        if rebuild_IV:
            self.build_IV()

    def build_IV(self, V=np.array([-0.1,0.1]), *args, **kwargs):
        self.IV_table = np.array([V, -self.IL*np.ones_like(V)])

    def __str__(self):
        return "Current Source: IL = " + self.get_value_text()
    
    def get_value_text(self):
        return f"{self.IL:.4f} A"
    def get_draw_func(self):
        return draw_CC_symbol

class Resistor(CircuitElement):
    def __init__(self, cond=1, tag=None):
        super().__init__(tag=tag)
        self.cond = cond
    def build_IV(self, V=np.array([-0.1,0.1]), *args, **kwargs):
        self.IV_table = np.array([V, V*self.cond])
    def set_cond(self,cond):
        self.cond = cond
        self.null_IV()
    def copy(self,source):
        self.set_cond(source.cond)
    def __str__(self):
        return "Resistor: R = " + self.get_value_text()
    def get_value_text(self):
        return f"{1/self.cond:.3f} ohm"
    def get_draw_func(self):
        return draw_resistor_symbol

class Diode(CircuitElement):
    def __init__(self,I0=1e-15,n=1,V_shift=0,tag=None,temperature=25): #V_shift is to shift the starting voltage, e.g. to define breakdown
        super().__init__(tag=tag)
        self.I0 = I0
        self.n = n
        self.V_shift = V_shift
        self.VT = get_VT(temperature)
        self.refI0 = I0
        self.refT = temperature

    def set_I0(self,I0):
        self.I0 = I0
        self.null_IV()

    def copy(self,source):
        self.n = source.n
        self.V_shift = source.V_shift
        self.VT = source.VT
        self.refI0 = source.refI0
        self.refT = source.refT
        self.set_I0(source.I0)

    def changeTemperature(self,temperature,rebuild_IV=True):
        self.VT = get_VT(temperature)
        old_ni  = 9.15e19*((self.refT+273.15)/300)**2*np.exp(-6880/(self.refT+273.15))
        new_ni  = 9.15e19*((temperature+273.15)/300)**2*np.exp(-6880/(temperature+273.15))
        scale_factor = (new_ni/old_ni)**(2/self.n)
        self.set_I0(self.refI0*scale_factor)
        if rebuild_IV:
            self.build_IV()
        
    def build_IV(self, V=None, max_num_points=100, *args, **kwargs):
        if V is None:
            if max_num_points is None:
                max_num_points = 100
            # assume that 2.0 A/cm2 is max you'll need
            if self.I0==0:
                Voc = 0.78
            else:
                Voc = self.n*self.VT*np.log(2.0/self.I0)
            V = [self.V_shift-1.1,self.V_shift-1.0,self.V_shift]+list(self.V_shift + Voc*np.log(np.arange(1,max_num_points))/np.log(max_num_points-1))
            V = np.array(V)
        I = self.I0*(np.exp((V-self.V_shift)/(self.n*self.VT))-1)
        self.IV_table = np.array([V,I])
    
class ForwardDiode(Diode):
    def __init__(self,I0=1e-15,n=1,tag=None): #V_shift is to shift the starting voltage, e.g. to define breakdown
        super().__init__(I0, n, V_shift=0,tag=tag)
    def build_IV(self, V=None, max_num_points=100, *args, **kwargs):
        super().build_IV(V,max_num_points)
    def __str__(self):
        return "Forward Diode: I0 = " + str(self.I0) + "A, n = " + str(self.n)
    def get_value_text(self):
        return f"I0 = {self.I0:.3e}A\nn = {self.n:.2f}"
    def get_draw_func(self):
        return draw_forward_diode_symbol
    
class PhotonCouplingDiode(ForwardDiode):
    def get_draw_func(self):
        return draw_LED_diode_symbol

class ReverseDiode(Diode):
    def __init__(self,I0=1e-15,n=1, V_shift=0,tag=None): #V_shift is to shift the starting voltage, e.g. to define breakdown
        super().__init__(I0, n, V_shift, tag=tag)
    def build_IV(self, V=None, max_num_points=100, *args, **kwargs):
        super().build_IV(V,max_num_points)
        self.IV_table[1,:] += self.I0
        self.IV_table *= -1
        self.IV_table = self.IV_table[:,::-1]
    def __str__(self):
        return "Reverse Diode: I0 = " + str(self.I0) + "A, n = " + str(self.n) + ", breakdown V = " + str(self.V_shift)
    def get_value_text(self):
        return f"I0 = {self.I0:.3e}A\nn = {self.n:.2f}\nbreakdown V = {self.V_shift:.2f}"
    def get_draw_func(self):
        return draw_reverse_diode_symbol

class CircuitGroup():
    def __init__(self,subgroups,connection="series",name=None,location=None,
                 rotation=0,x_mirror=1,y_mirror=1,extent=None):
        self.connection = connection
        self.subgroups = subgroups
        for element in self.subgroups:
            element.parent = self
        self.parent = None
        self.IV_table = None
        self.dark_IV_table = None
        self.name = name
        if location is None:
            self.location = np.array([0,0])
        else:
            self.location = location
        self.rotation = rotation
        self.x_mirror = x_mirror
        self.y_mirror = y_mirror
        if extent is not None:
            self.extent = extent
        else:
            self.extent = get_extent(subgroups)
        self.circuit_diagram_extent = get_circuit_diagram_extent(subgroups,connection)
        self.operating_point = None #V,I
        self.aux = {}
    
    def null_IV(self, keep_dark=False):
        if self.IV_table is not None or self.dark_IV_table is not None:
            self.IV_table = None
            if keep_dark==False:
                self.dark_IV_table = None
            if self.parent is not None:
                self.parent.null_IV(keep_dark=keep_dark)

    def null_all_IV(self):
        self.IV_table = None
        if hasattr(self,"dark_IV_table"):
            self.dark_IV_table = None
        for element in self.subgroups:
            if isinstance(element,CircuitElement):
                element.IV_table = None
            else:
                element.null_all_IV()

    def reassign_parents(self):
        for element in self.subgroups:
            element.parent = self
            if isinstance(element,CircuitGroup):
                element.reassign_parents()

    def set_operating_point(self,V=None,I=None):
        if V is not None:
            I = np.interp(V,self.IV_table[0,:],self.IV_table[1,:])
        elif I is not None:
            V = np.interp(I,self.IV_table[1,:],self.IV_table[0,:])
        for element in self.subgroups:
            if self.connection == "series": # then all elements have same current
                target_I = I
                # solar cell needs to scale IV table by area
                if hasattr(self,"shape") and self.area is not None:
                    target_I /= self.area
                element.set_operating_point(V=None,I=target_I)
            else: # then all elements have same voltage
                element.set_operating_point(V=V,I=None)
        self.operating_point = [V,I]
        # cells also store Vint
        if hasattr(self,"shape"):
            self.operating_point.append(self.diode_branch.operating_point[0])

    def removeElementOfTag(self,tag):
        for element in self.subgroups[:]:
            if isinstance(element,CircuitElement):
                if element.tag==tag:
                    self.subgroups.remove(element)
            elif isinstance(element,CircuitGroup):
                element.removeElementOfTag(tag)
        self.null_IV()

    def set_temperature(self,temperature,rebuild_IV=True):
        diodes = self.findElementType(Diode)
        for diode in diodes:
            diode.changeTemperature(temperature,rebuild_IV=False)
        currentSources = self.findElementType(CurrentSource)
        for currentSource in currentSources:
            currentSource.changeTemperatureAndSuns(temperature=temperature,rebuild_IV=False)
        if rebuild_IV:
            self.build_IV()

    def findElementType(self,type,serialize=False,path=[]):
        list_ = []
        for i, element in enumerate(self.subgroups):
            if isinstance(element,type):
                list_.append(element)
            elif isinstance(element,CircuitGroup):
                list_.extend(element.findElementType(type,serialize=serialize))
        if serialize:
            for i, element in enumerate(list_):
                element.name = str(i)
        return list_
    
    def build_IV(self, max_num_points=None, cap_current=None):
        # if solar cell, then express in current density
        Vints = None
        if hasattr(self,"shape") and self.area is not None and cap_current is not None:
            cap_current /= self.area
        all_circuit_element_children = True
        for element in self.subgroups:
            if isinstance(element,CircuitGroup):
                all_circuit_element_children = False
            if element.IV_table is None:
                element.build_IV(max_num_points=max_num_points,cap_current=cap_current)
        shift_IV_only = False 
        total_IL = 0.0       
        if self.connection=="series":
            # add voltage
            Is = []
            for element in self.subgroups:
                Is.extend(list(element.IV_table[1,:]))
            Is = np.sort(np.array(Is))
            Is = np.unique(Is)
            Vs = np.zeros_like(Is)
            Vints = np.zeros_like(Is)
            # do reverse order to allow for photon coupling
            pc_IVs = []
            for element in reversed(self.subgroups):
                IV_table = element.IV_table
                for pc_IV in pc_IVs:
                    IV_table[1,:] += interp_(IV_table[0,:], pc_IV[0,:], pc_IV[1,:])
                V = interp_(Is,IV_table[1,:],IV_table[0,:])
                Vs += V
                if hasattr(self,"shape") and isinstance(element,CircuitGroup):
                    Vints += V
                pc_IVs = []
                if hasattr(element,"photon_coupling_diodes"):
                    for pc in element.photon_coupling_diodes:
                        pc_IVs.append(pc.IV_table)
                if element.IV_table.shape[0]==3 and np.max(element.IV_table[2,:])>0:
                    Vint = interp_(Is,IV_table[2,:],IV_table[0,:])
                    Vints += Vint
        else:
            # add current
            for element in self.subgroups:
                if isinstance(element,CurrentSource):
                    total_IL -= element.IL
            if self.dark_IV_table is not None:
                shift_IV_only = True
                self.IV_table = self.dark_IV_table.copy()
                if hasattr(self,"shape") and self.area is not None:
                    total_IL *= self.area
                self.IV_table[1,:] += total_IL
            else:
                Vs = []
                for element in self.subgroups:
                    Vs.extend(list(element.IV_table[0,:]))
                Vs = np.sort(np.array(Vs))
                Vs = np.unique(Vs)
                Is = np.zeros_like(Vs)
                for element in self.subgroups:
                    if not isinstance(element,CurrentSource):
                        Is += interp_(Vs,element.IV_table[0,:],element.IV_table[1,:])
                Is += total_IL
                if cap_current is not None:
                    find_ = np.where(np.abs(Is) < cap_current)[0]
                    Vs = Vs[find_]
                    Is = Is[find_]

        if shift_IV_only==False:
            self.IV_table = np.array([Vs,Is])
            if max_num_points is None:
                self.IV_table = np.array([Vs,Is])
            else:
                V_range = np.max(Vs)-np.min(Vs)
                I_range = np.max(Is)-np.min(Is)
                V_segments = Vs[1:]-Vs[:-1]
                I_segments = Is[1:]-Is[:-1]
                segment_lengths = np.sqrt((V_segments/V_range)**2+(I_segments/I_range)**2)
                total_length = np.sum(segment_lengths)
                ideal_segment_length = total_length / max_num_points
                short_segments = np.where(segment_lengths < ideal_segment_length)[0]
                long_segments = np.where(segment_lengths >= ideal_segment_length)[0]
                if len(short_segments)>0:
                    short_segment_lengths = segment_lengths[short_segments]
                    short_segment_lengths_cum = np.cumsum(short_segment_lengths)
                    short_segment_left_V = Vs[short_segments]
                    long_segment_left_V = Vs[long_segments]
                    ideal_Vs = np.linspace(0,short_segment_lengths_cum[-1],max_num_points-len(long_segments))
                    short_segment_lengths_cum -= short_segment_lengths_cum[0]
                    index = np.searchsorted(short_segment_lengths_cum, ideal_Vs, side='right')-1
                    new_Vs = short_segment_left_V[index] + ideal_Vs - short_segment_lengths_cum[index]
                    new_Vs = list(new_Vs)
                    new_Vs.extend(long_segment_left_V)
                    new_Vs = np.sort(np.array(new_Vs))
                    if new_Vs[0] > Vs[0]:
                        new_Vs = np.hstack((Vs[0],new_Vs))
                    if new_Vs[-1] < Vs[-1]:
                        new_Vs = np.hstack((new_Vs,Vs[-1]))
                    new_Is = interp_(new_Vs,Vs,Is)
                    if Vints is not None:
                        new_Vints = interp_(new_Is,Is,Vints)
                        self.IV_table = np.array([new_Vs, new_Is, new_Vints])
                    else:
                        self.IV_table = np.array([new_Vs, new_Is])
                else:
                    if Vints is not None:
                        self.IV_table = np.array([Vs,Is,Vints])
                    else:
                        self.IV_table = np.array([Vs,Is])
            if all_circuit_element_children:
                self.dark_IV_table = self.IV_table.copy()
                self.dark_IV_table[1,:] -= total_IL
            # solar cell needs to scale IV table by area
            if hasattr(self,"shape") and self.area is not None:
                self.IV_table[1,:] *= self.area
                if self.dark_IV_table is not None:
                    self.dark_IV_table[1,:] *= self.area
    
    def __str__(self):
        word = self.connection + " connection:\n"
        for i, element in enumerate(self.subgroups):
            if isinstance(element,CircuitGroup):
                word += "Subgroup " + str(i) + ":\n"
            word += str(element) + "\n"
        return word    
    
    def draw(self, ax=None, x=0, y=0, display_value=False):
        global pbar
        draw_immediately = False
        if ax is None:
            num_of_elements = len(self.findElementType(CircuitElement))
            pbar = tqdm(total=num_of_elements)
            fig, ax = plt.subplots()
            draw_immediately = True
        
        current_x = x - self.circuit_diagram_extent[0]/2
        current_y = y - self.circuit_diagram_extent[1]/2
        if self.connection != "series":
            current_y += 0.1
        for i, element in enumerate(self.subgroups):
            if isinstance(element,CircuitElement):
                pbar.update(1)
            center_x = current_x+element.circuit_diagram_extent[0]/2
            center_y = current_y+element.circuit_diagram_extent[1]/2
            if self.connection == "series":
                center_x = x
            else:
                center_y = y
            element.draw(ax=ax, x=center_x, y=center_y, display_value=display_value)
            if self.connection=="series":
                if i > 0:
                    line = plt.Line2D([x,x],[current_y-y_spacing, current_y], color="black", linewidth=1.5)
                    ax.add_line(line)
                current_y += element.circuit_diagram_extent[1]+y_spacing
            else:
                line = plt.Line2D([center_x,center_x], [center_y+element.circuit_diagram_extent[1]/2,y+self.circuit_diagram_extent[1]/2], color="black", linewidth=1.5)
                ax.add_line(line)
                line = plt.Line2D([center_x,center_x], [center_y-element.circuit_diagram_extent[1]/2,y-self.circuit_diagram_extent[1]/2], color="black", linewidth=1.5)
                ax.add_line(line)
                if i > 0:
                    line = plt.Line2D([center_x,current_x-x_spacing-self.subgroups[i-1].circuit_diagram_extent[0]/2], [y+self.circuit_diagram_extent[1]/2,y+self.circuit_diagram_extent[1]/2], color="black", linewidth=1.5)
                    ax.add_line(line)
                    line = plt.Line2D([center_x,current_x-x_spacing-self.subgroups[i-1].circuit_diagram_extent[0]/2], [y-self.circuit_diagram_extent[1]/2,y-self.circuit_diagram_extent[1]/2], color="black", linewidth=1.5)
                    ax.add_line(line)
                current_x += element.circuit_diagram_extent[0]+x_spacing
        if draw_immediately:
            pbar.close()
            line = plt.Line2D([x,x], [y-self.circuit_diagram_extent[1]/2,y-self.circuit_diagram_extent[1]/2-0.2], color="black", linewidth=1.5)
            ax.add_line(line)
            line = plt.Line2D([x,x], [y+self.circuit_diagram_extent[1]/2,y+self.circuit_diagram_extent[1]/2+0.2], color="black", linewidth=1.5)
            ax.add_line(line)
            draw_symbol(draw_earth_symbol, ax=ax,  x=x, y=y-self.circuit_diagram_extent[1]/2-0.3)
            draw_symbol(draw_pos_terminal_symbol, ax=ax,  x=x, y=y+self.circuit_diagram_extent[1]/2+0.25)
            ax.set_aspect('equal')
            ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
            for spine in ax.spines.values():
                spine.set_visible(False)
            fig.tight_layout()
            plt.show()

def get_extent(elements, center=True):
    x_bounds = [None,None]
    y_bounds = [None,None]
    for element in elements:
        if hasattr(element,"extent") and element.extent is not None:
            xs = [element.location[0]-element.extent[0]/2,element.location[0]+element.extent[0]/2]
            ys = [element.location[1]-element.extent[1]/2,element.location[1]+element.extent[1]/2]
            if x_bounds[0] is None:
                x_bounds[0] = xs[0]
            else:
                x_bounds[0] = min(x_bounds[0],xs[0])
            if x_bounds[1] is None:
                x_bounds[1] = xs[1]
            else:
                x_bounds[1] = max(x_bounds[1],xs[1])
            if y_bounds[0] is None:
                y_bounds[0] = ys[0]
            else:
                y_bounds[0] = min(y_bounds[0],ys[0])
            if y_bounds[1] is None:
                y_bounds[1] = ys[1]
            else:
                y_bounds[1] = max(y_bounds[1],ys[1])
    if (x_bounds[0] is not None) and (x_bounds[1] is not None) and (y_bounds[0] is not None) and (y_bounds[1] is not None):
        if center:
            center = [0.5*(x_bounds[0]+x_bounds[1]),0.5*(y_bounds[0]+y_bounds[1])]
            for element in elements:
                if hasattr(element,"extent") and element.extent is not None:
                    element.location[0] -= center[0]
                    element.location[1] -= center[1]
        return [x_bounds[1]-x_bounds[0],y_bounds[1]-y_bounds[0]]
    else:
        return None

def get_circuit_diagram_extent(elements,connection):
    total_extent = [0.0,0.0]
    for i, element in enumerate(elements):
        extent_ = element.circuit_diagram_extent
        if connection=="series":
            total_extent[0] = max(total_extent[0], extent_[0])
            total_extent[1] += extent_[1]
            if i > 0:
                total_extent[1] += y_spacing
        else:
            total_extent[1] = max(total_extent[1], extent_[1])
            total_extent[0] += extent_[0]
            if i > 0:
                total_extent[0] += x_spacing
    if connection!="series":
        total_extent[1] += 0.2 # the connectors
    return total_extent

def tile_elements(elements, rows=None, cols=None, x_gap = 0.0, y_gap = 0.0, turn=True):
    assert((rows is not None) or (cols is not None))
    if rows is None:
        rows = int(np.ceil(float(len(elements))/float(cols)))
    if cols is None:
        cols = int(np.ceil(float(len(elements))/float(rows)))
    row = 0
    col = 0
    rotation = 0
    pos = np.array([0,0]).astype(float)
    bounds = np.array([0,0]).astype(float)
    max_x_extent = 0.0
    for element in elements:
        x_extent = element.extent[0]
        max_x_extent = max(max_x_extent,x_extent)
        y_extent = element.extent[1]
        element.location = pos.copy()
        element.rotation = rotation
        row += 1
        if col == cols-1:
            bounds[0] = max(bounds[0],pos[0] + max_x_extent)
        if row == rows:
            bounds[1] = max(bounds[1],pos[1] + y_extent)
        if row < rows:
            if rotation==0:
                pos[1] += y_extent + y_gap
            else:
                pos[1] -= (y_extent + y_gap)
        else:
            row = 0
            col += 1
            pos[0] += max_x_extent + x_gap
            max_x_extent = 0.0
            if turn:
                rotation = 180 - rotation
            else:
                pos[1] = 0

def circuit_deepcopy(circuit_group):
    circuit_group2 = copy.deepcopy(circuit_group)
    circuit_group2.reassign_parents()
    return circuit_group2
