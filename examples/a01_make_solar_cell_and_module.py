from PV_Circuit_Model.cell import *
from PV_Circuit_Model.module import *
from PV_Circuit_Model.cell_analysis import *

def main(display=True):
    # desired cell parameters
    width = 18.2 # cm
    Jsc = 0.042 # A/cm2
    Voc = 0.735 # V
    FF = 0.82 
    breakdown_V = -10 # V
    J0_rev = 100e-15 # A/cm2
    Rshunt = 10000 # ohm-cm2
    Rs = 0.3333 # ohm-cm2
    thickness = 150e-4 #um

    # butterfly module layout
    n_cells = [22,6]
    num_cells_per_halfstring = n_cells[0]
    num_half_strings = n_cells[1]

    # make the solar cell 
    shape, area = wafer_shape(width/2, width)
    J01, J02 = estimate_cell_J01_J02(Jsc,Voc,FF=FF,Rs=Rs,Rshunt=Rshunt,thickness=thickness)
    cell = make_solar_cell(Jsc, J01, J02, Rshunt, Rs, area, shape, thickness, breakdown_V, J0_rev)

    if display:
        # draw its circuit model representation
        cell.draw(display_value=True)
        # plot its IV curve
        cell.plot()
        cell.show()
        # write out its constituent parts and values
        print(cell)

    cells = [circuit_deepcopy(cell) for _ in range(num_half_strings*num_cells_per_halfstring)]
    module = make_butterly_module(cells, num_strings=num_half_strings // 2, num_cells_per_halfstring=num_cells_per_halfstring)
    if display:
        # draw module cells layout
        module.draw_cells(show_names=True)
        # draw its circuit model representation
        module.draw()
        # plot its IV curve
        module.plot()
        module.show()
        # write out its constituent parts and values
        print(module)
        
    return module

if __name__ == "__main__": 
    main()
