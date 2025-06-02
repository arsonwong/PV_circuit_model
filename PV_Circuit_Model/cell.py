import numpy as np
from PV_Circuit_Model.circuit_model import *
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# For Si, the intrinsic recomb limit for every cm3, for 1e15cm-3 doped n type,
# works out to be around
# I01 = 9.664293091887626e-14A, n1 = 1
# I02 = 4.882772711642952e-19A, n2 = 0.6849
class instrinsic_Si():
    J01 = 9.664293091887626e-14
    n1 = 1
    J0x = 4.882772711642952e-19
    nx = 0.6849
    Jsc_fractional_temp_coeff = 0.0004

class Cell(CircuitGroup):
    def __init__(self,components,connection="series",area=None,location=None,
                 rotation=0,shape=None,thickness=None,name=None,temperature=25,Suns=1.0):
        x_extent = 0.0
        y_extent = 0.0
        if shape is not None:
            x_extent = np.max(shape[:,0])-np.min(shape[:,0])
            y_extent = np.max(shape[:,1])-np.min(shape[:,1])
        super().__init__(components, connection,location=location,rotation=rotation,
                         name=name,extent=np.array([x_extent,y_extent]).astype(float))
        self.area = area
        self.shape = shape
        self.build_IV()
        self.temperature = temperature
        self.set_temperature(temperature)
        self.Suns = Suns
        self.thickness = thickness
        self.get_branches()
        self.photon_coupling_diodes = self.findElementType(PhotonCouplingDiode)

    def get_branches(self):
        if self.connection=="series":
            for branch in self.subgroups:
                if isinstance(branch,Resistor):
                    self.series_resistor = branch
                else:
                    self.diode_branch = branch
        else:
            self.series_resistor = None
            self.diode_branch = self

    # a weak copy, only the parameters
    def copy(self,cell2):
        self.thickness = cell2.thickness
        self.temperature = cell2.temperature
        self.Suns = cell2.Suns
        if self.series_resistor is not None:
            self.series_resistor.copy(cell2.series_resistor)
        for i, element in enumerate(self.diode_branch.subgroups):
            if i < len(cell2.diode_branch.subgroups):
                element.copy(cell2.diode_branch.subgroups[i])

    def set_Suns(self,Suns,rebuild_IV=True):
        self.Suns = Suns
        currentSources = self.findElementType(CurrentSource)
        for currentSource in currentSources:
            currentSource.changeTemperatureAndSuns(Suns=Suns)
        if rebuild_IV:
            self.build_IV()
    
    def set_temperature(self,temperature,rebuild_IV=True):
        super().set_temperature(temperature)
        self.temperature = temperature
        if rebuild_IV:
            self.build_IV()

    def JL(self):
        JL = 0.0
        currentSources = self.findElementType(CurrentSource)
        for currentSource in currentSources:
            JL += currentSource.IL
        return JL

    def IL(self):
        return self.JL()*self.area     
    
    def set_JL(self,JL,Suns=1.0,temperature=25,rebuild_IV=True):
        currentSources = self.findElementType(CurrentSource)
        for currentSource in currentSources:
            if currentSource.tag != "defect":
                currentSource.refSuns = Suns
                currentSource.refIL = JL
                currentSource.refT = temperature
                currentSource.changeTemperatureAndSuns(
                    temperature=self.temperature,Suns=self.Suns,rebuild_IV=False)
                break
        if rebuild_IV:
            self.build_IV()

    def set_IL(self,IL,Suns=1.0,temperature=25,rebuild_IV=True):
        self.set_JL(IL/self.area,Suns=Suns,temperature=temperature,rebuild_IV=rebuild_IV)
    
    def J0(self,n):
        J0 = 0.0
        diodes = self.findElementType(ForwardDiode)
        for diode in diodes:
            if diode.n==n:
                if diode.tag != "intrinsic" and not isinstance(diode,PhotonCouplingDiode):
                    J0 += diode.I0
        return J0
    def J01(self):
        return self.J0(n=1)
    def J02(self):
        return self.J0(n=2)     
    
    def PC_J0(self,n):
        J0 = 0.0
        diodes = self.findElementType(PhotonCouplingDiode)
        for diode in diodes:
            if diode.n==n:
                J0 += diode.I0
        return J0
    def PC_J01(self):
        return self.PC_J0(n=1)
    def PC_I0(self,n):
        return self.PC_J0(n)*self.area
    def PC_I01(self):
        return self.PC_I0(n=1)
    
    def I0(self,n):
        return self.J0(n)*self.area
    def I01(self):
        return self.I0(n=1)
    def I02(self):
        return self.I0(n=2)    
    
    def set_J0(self,J0,n,temperature=25,rebuild_IV=True):
        diodes = self.findElementType(ForwardDiode)
        for diode in diodes:
            if diode.tag != "defect" and diode.tag != "intrinsic" and diode.n==n and not isinstance(diode,PhotonCouplingDiode):
                diode.refI0 = J0
                diode.refT = temperature
                diode.changeTemperature(temperature=self.temperature,rebuild_IV=False)
                break
        if rebuild_IV:
            self.build_IV()
    def set_J01(self,J0,temperature=25,rebuild_IV=True):
        self.set_J0(J0,n=1,temperature=temperature,rebuild_IV=rebuild_IV)
    def set_J02(self,J0,temperature=25,rebuild_IV=True):
        self.set_J0(J0,n=2,temperature=temperature,rebuild_IV=rebuild_IV)

    def set_I0(self,I0,n,temperature=25,rebuild_IV=True):
        self.set_J0(I0/self.area, n=n, temperature=temperature,rebuild_IV=rebuild_IV)
    def set_I01(self,I0,temperature=25,rebuild_IV=True):
        self.set_I0(I0,n=1,temperature=temperature,rebuild_IV=rebuild_IV)
    def set_I02(self,I0,temperature=25,rebuild_IV=True):
        self.set_I0(I0,n=2,temperature=temperature,rebuild_IV=rebuild_IV)

    def set_PC_J0(self,J0,n,temperature=25,rebuild_IV=True):
        diodes = self.findElementType(PhotonCouplingDiode)
        for diode in diodes:
            if diode.tag != "defect" and diode.tag != "intrinsic" and diode.n==n:
                diode.refI0 = J0
                diode.refT = temperature
                diode.changeTemperature(temperature=self.temperature,rebuild_IV=False)
                break
        if rebuild_IV:
            self.build_IV()
    def set_PC_J01(self,J0,temperature=25,rebuild_IV=True):
        self.set_PC_J0(J0,n=1,temperature=temperature,rebuild_IV=rebuild_IV)
    def set_PC_I0(self,I0,n,temperature=25,rebuild_IV=True):
        self.set_PC_J0(I0/self.area, n=n, temperature=temperature,rebuild_IV=rebuild_IV)
    def set_PC_I01(self,I0,temperature=25,rebuild_IV=True):
        self.set_PC_I0(I0,n=1,temperature=temperature,rebuild_IV=rebuild_IV)
    
    def specific_Rs_cond(self):
        if self.series_resistor is None:
            return np.inf
        return self.series_resistor.cond
    def Rs_cond(self):
        return self.specific_Rs_cond()*self.area
    def specific_Rs(self):
        return 1/self.specific_Rs_cond()
    def Rs(self):
        return 1/self.Rs_cond()
    
    def set_specific_Rs_cond(self,cond):
        if self.series_resistor is not None:
            self.series_resistor.set_cond(cond)
    def set_Rs_cond(self,cond):
        self.set_specific_Rs_cond(cond/self.area)
    def set_specific_Rs(self,Rs):
        self.set_specific_Rs_cond(1/Rs)
    def set_Rs(self,Rs):
        self.set_specific_Rs(Rs*self.area)
    
    def specific_shunt_cond(self):
        Rsh_cond = 0.0
        shunt_resistors = self.diode_branch.findElementType(Resistor)
        for res in shunt_resistors:
            Rsh_cond += res.cond
        return Rsh_cond
    def shunt_cond(self):
        return self.specific_shunt_cond()*self.area
    def specific_shunt_res(self):
        return 1/self.specific_shunt_cond()
    def shunt_res(self):
        return 1/self.shunt_cond()
    
    def set_specific_shunt_cond(self,cond):
        shunt_resistors = self.diode_branch.findElementType(Resistor)
        for res in shunt_resistors:
            if res.tag != "defect":
                res.set_cond(cond)
                break
    def set_shunt_cond(self,cond):
        self.set_specific_shunt_cond(cond/self.area)
    def set_specific_shunt_res(self,Rsh):
        self.set_specific_shunt_cond(1/Rsh)
    def set_shunt_res(self,Rsh):
        self.set_specific_shunt_res(Rsh*self.area)

class Module(CircuitGroup):
    def __init__(self,subgroups,connection="series",location=np.array([0,0]),
                 rotation=0,cap_current=None,name=None,temperature=25,Suns=1.0):
        super().__init__(subgroups, connection,location=location,rotation=rotation,name=name)
        if self.location is None:
            self.location = np.array([0,0])
        self.cap_current = cap_current
        cells = self.findElementType(Cell,serialize=True)
        self.cells = cells
        self.temperature = temperature
        self.set_temperature(temperature)
        self.Suns = Suns
        self.set_Suns(Suns)     
    def set_Suns(self,Suns, rebuild_IV=True):
        for cell in self.cells:
            cell.set_Suns(Suns=Suns, rebuild_IV=False)
        if rebuild_IV:
            self.build_IV()
    def set_temperature(self,temperature, rebuild_IV=True):
        super().set_temperature(temperature,rebuild_IV=False)
        self.temperature = temperature
        if rebuild_IV:
            self.build_IV()
    def build_IV(self, max_num_points=500):
        super().build_IV(max_num_points=max_num_points,
                         cap_current=self.cap_current)
    
def draw_cells(self: CircuitGroup,display=True,show_names=False,colour_what="Vint"):
    shapes = []
    names = []
    Vints = []
    EL_Vints = []
    if hasattr(self,"shape"): # a solar cell
        shapes.append(self.shape.copy())
        names.append(self.name)
        if self.operating_point is not None and len(self.operating_point)==3:
            Vints.append(self.operating_point[2])
        if self.aux is not None and "EL_Vint" in self.aux:
            EL_Vints.append(self.aux["EL_Vint"])
    else:
        for element in self.subgroups:
            if hasattr(element,"extent") and element.extent is not None:
                shapes_, names_, Vints_, EL_Vints_ = element.draw_cells(display=False)
                shapes.extend(shapes_)
                names.extend(names_)
                Vints.extend(Vints_)
                EL_Vints.extend(EL_Vints_)

    has_Vint = False
    has_EL_Vint = False
    if len(EL_Vints)==len(shapes) and colour_what=="EL_Vint": # every cell has a EL_Vint
        has_EL_Vint = True
        norm = mcolors.Normalize(vmin=min(EL_Vints), vmax=max(EL_Vints))
        cmap = cm.viridis
    elif len(Vints)==len(shapes) and colour_what=="Vint": # every cell has a Vint
        has_Vint = True
        norm = mcolors.Normalize(vmin=min(Vints), vmax=max(Vints))
        cmap = cm.viridis 
            
    for i, shape in enumerate(shapes):
        cos = np.cos(np.pi/180*self.rotation)
        sin = np.sin(np.pi/180*self.rotation)
        new_shape = shape.copy()
        new_shape[:,0] = shape[:,0]*cos + shape[:,1]*sin
        new_shape[:,1] = shape[:,1]*cos - shape[:,0]*sin
        if self.x_mirror == -1:
            new_shape[:,0] *= -1
        if self.y_mirror == -1:
            new_shape[:,1] *= -1
        new_shape[:,0] += self.location[0]
        new_shape[:,1] += self.location[1]

        shapes[i] = new_shape
    if display:
        fig, ax = plt.subplots()
        for i, shape in enumerate(shapes):
            color = 'skyblue'
            if has_EL_Vint:
                color = cmap(norm(EL_Vints[i]))
            elif has_Vint:
                color = cmap(norm(Vints[i]))
            polygon = patches.Polygon(shape, closed=True, facecolor=color, edgecolor='black')
            x = 0.5*(np.max(shape[:,0])+np.min(shape[:,0]))
            y = 0.5*(np.max(shape[:,1])+np.min(shape[:,1]))
            if show_names:
                ax.text(x, y, names[i], fontsize=8, color='black')
            ax.add_patch(polygon)
        ax.set_xlim(self.location[0]-self.extent[0]/2*1.1, self.location[0]+self.extent[0]/2*1.1)
        ax.set_ylim(self.location[1]-self.extent[1]/2*1.1, self.location[1]+self.extent[1]/2*1.1)
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        fig.tight_layout()
        ax.set_aspect('equal')
        plt.show()
    return shapes, names, Vints, EL_Vints
CircuitGroup.draw_cells = draw_cells

def wafer_shape(L, W, ingot_center=None, ingot_diameter=None):
    shape = np.array([[-W/2,-L/2],[W/2,-L/2],[W/2,L/2],[-W/2,L/2]])
    x = shape[:,0]
    y = shape[:,1]
    area = 0.5 * np.abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    return shape, area

# note: always made at 25C 1 Sun
def make_solar_cell(Jsc=0.042, J01=10e-15, J02=2e-9, Rshunt=1e6, Rs=0.0, area=1.0, 
                    shape=None, thickness=180e-4, breakdown_V=-10, J0_rev=100e-15,
                    J01_photon_coupling=0.0, Si_intrinsic_limit=True):
    elements = [CurrentSource(IL=Jsc, temp_coeff = instrinsic_Si.Jsc_fractional_temp_coeff*Jsc),
                ForwardDiode(I0=J01,n=1),
                ForwardDiode(I0=J02,n=2)]
    if J01_photon_coupling > 0:
        elements.append(PhotonCouplingDiode(I0=J01_photon_coupling,n=1))
    if Si_intrinsic_limit:
        elements.extend([ForwardDiode(I0=instrinsic_Si.J01*thickness, n=instrinsic_Si.n1,tag="intrinsic"),
                ForwardDiode(I0=instrinsic_Si.J0x*thickness, n=instrinsic_Si.nx,tag="intrinsic")])
    elements.extend([ReverseDiode(I0=J0_rev, n=1, V_shift = -breakdown_V),
                Resistor(cond=1/Rshunt)])
    if Rs == 0.0:
        cell = Cell(elements,"parallel",area=area,thickness=thickness,location=np.array([0.0,0.0]).astype(float),shape=shape,name="cell")
    else:
        group = CircuitGroup(elements,"parallel")
        cell = Cell([group,Resistor(cond=1/Rs)],"series",area=area,thickness=thickness,location=np.array([0.0,0.0]).astype(float),shape=shape,name="cell")
    
    return cell
