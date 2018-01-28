
from my_models import mason_model
import OpenPNM
import numpy as np
import OpenPNM as op
from my_custom_phases import Brine
from my_custom_phases import C02
import new_models
import xlwt
from xlrd import open_workbook

from xlutils.copy import copy
rock=['berea-3.20','bentheimer-3.18','mtsimon-2.55','mtsimon-2.80']

extracted_network = 'berea-3.20'
pn = op.Utilities.IO.Statoil.load(path=r'./', prefix=extracted_network)

pn.name = extracted_network
op.export_data(network=pn, filename=extracted_network)
pn.trim(pn.check_network_health()['trim_pores'])
print(op)
"""
with cls._read_file(filename=extracted_network+"_node1", ext='dat') as f:
                row_0 = f.readline().split()
                num_lines = int(row_0[0])

"""
## define geom and put ALL the data which stored on the Network object into geom
geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts,name='internal')
for item in pn.props():
   if item not in ['throat.conns', 'pore.coords']:
       #geom.update({item: pn.pop(item)})
        geom.update({item: pn[item][~np.isnan(pn[item])]})
#print(geom)


geom['pore.diameter'] = 2*geom['pore.radius']
geom.models.add(propname='pore.area',
                model=op.Geometry.models.pore_area.spherical)
geom['throat.diameter'] = 2*geom['throat.radius']

geom.models.add(propname='throat.area',
                model=op.Geometry.models.throat_area.cylinder)


#print(pn,"after adding boundary")


Wphase = Brine(network=pn)
Nphase = C02(network=pn)
Nphase['throat.contact_angle'] = 30
Nphase['throat.surface_tension'] = 0.072

Wphase['throat.contact_angle'] = 150
Wphase['throat.surface_tension'] = 0.072



print(Wphase,"Wphase")
print(Nphase,"Nphase")

phys_Nphase_internal = op.Physics.GenericPhysics(network=pn, phase=Nphase, geometry=geom)
phys_Wphase_internal = op.Physics.GenericPhysics(network=pn, phase=Wphase, geometry=geom)
#phys_Nphase_boundary = op.Physics.GenericPhysics(network=pn, phase=Nphase, geometry=boun)
#phys_Wphase_boundary = op.Physics.GenericPhysics(network=pn, phase=Wphase, geometry=boun)

phys_Nphase_internal.models.add(propname='throat.capillary_pressure',model=mason_model)

phys_Nphase_internal.models.add(propname='throat.hydraulic_conductance',model=new_models.noncircular_capillary)


phys_Wphase_internal.models.add(propname='throat.capillary_pressure',model=mason_model)

phys_Wphase_internal.models.add(propname='throat.hydraulic_conductance',model=new_models.noncircular_capillary)


"""
sorted(list(pn.geometries.keys()))
temp = [item.regenerate for item in pn.physics.values()]
"""
#print(phys_Nphase_internal.models['throat.capillary_pressure']['f'])  # Inspect present value

#phys_Nphase_internal.models['throat.capillary_pressure']['f'] = 0.75  # Change value
#phys_Nphase_internal.models.regenerate()  # Regenerate model with new 'f' value

## the unit of capillary pressure is dependent on the Washburn's equation, surface tension(N/m)  / radius
inv = op.Algorithms.Drainage(network=pn,)
inv.setup(invading_phase=Nphase, defending_phase=Wphase)
inv.set_inlets(pores=pn.pores('*lets'))
inv.run(npts=150)

MIPdata = inv.get_drainage_data()
MIParray = np.array([MIPdata[key] for key in ('capillary_pressure', 'invading_phase_saturation', 'defending_phase_saturation')])
Pc = MIParray[0][:]
Snw=MIParray[1][:]
Sw = MIParray[2][:]


T_number_N =len(pn["throat.all"])
T_inv_N=np.arange(1,T_number_N+1,1)
P_number_N = len(pn["pore.all"])
P_inv_N=np.arange(1,P_number_N+1,1)
## calculate K_relative for increasing pressure
K_rel_array=[]
for item in Pc:
    #print(item)
   # Ts = pn.throats('boundary', mode='not')
    Pi = inv['pore.inv_Pc'] <= item
   # Ti = inv['throat.inv_Pc'][Ts] <= item
    Ti = inv['throat.inv_Pc'] <= item

    phys_Nphase_internal['throat.hydraulic_conductance'][~Ti] = 1e-20

    Nphase_flow = op.Algorithms.StokesFlow(network=pn, phase=Nphase)
    Nphase_flow.set_boundary_conditions(pores=pn.pores('inlets'), bcvalue=200000, bctype='Dirichlet')
    Nphase_flow.set_boundary_conditions(pores=pn.pores('outlets'), bcvalue=100000, bctype='Dirichlet')
    Nphase_flow.run()
    Q_partial = Nphase_flow.rate(pores=pn.pores('outlets'))

## find the absolute permeability of water
    phys_Nphase_internal.models.add(propname='throat.hydraulic_conductance',model=new_models.noncircular_capillary)
    phys_Nphase_internal.models.regenerate()
    Nphase_flow.run()
    Q_full = Nphase_flow.rate(pores=pn.pores('outlets'))
    K_rel = Q_partial/Q_full
    K_rel_array.extend(K_rel)

    T_inv_N=np.vstack([T_inv_N,~Ti])
    P_inv_N=np.vstack([P_inv_N,~Pi])
np.savetxt('T_inv_N.txt', T_inv_N.transpose(),fmt="%5i")
np.savetxt('p_inv_N.txt', P_inv_N.transpose(),fmt="%5i")

###################################################################
T_number_W =len(pn["throat.all"])
T_inv_W=np.arange(1,T_number_W+1,1)
P_number_W = len(pn["pore.all"])
P_inv_W=np.arange(1,P_number_W+1,1)


K_rel_array1=[]
for item in Pc:
    Pi = inv['pore.inv_Pc'] <= item
    Ti = inv['throat.inv_Pc'] <= item
    phys_Wphase_internal['throat.hydraulic_conductance'][Ti] = 1e-20
    Wphase_flow = op.Algorithms.StokesFlow(network=pn, phase=Wphase)
    Wphase_flow.set_boundary_conditions(pores=pn.pores('inlets'), bcvalue=200000, bctype='Dirichlet')
    Wphase_flow.set_boundary_conditions(pores=pn.pores('outlets'), bcvalue=100000, bctype='Dirichlet')
    Wphase_flow.run()
    Q_partial = Wphase_flow.rate(pores=pn.pores('outlets'))

    phys_Wphase_internal.models.add(propname='throat.hydraulic_conductance',model=new_models.noncircular_capillary)
    phys_Wphase_internal.models.regenerate()
    Wphase_flow.run()
    Q_full = Wphase_flow.rate(pores=pn.pores('outlets'))
    K_rel = Q_partial/Q_full
    K_rel_array1.extend(K_rel)

    T_inv_W=np.vstack([T_inv_W,Ti])
    P_inv_W=np.vstack([P_inv_W,Pi])
np.savetxt('T_inv_W.txt', T_inv_W.transpose(),fmt="%5i")
np.savetxt('p_inv_W.txt', P_inv_W.transpose(),fmt="%5i")




"""
plt.figure(1, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')
plt.plot(Sw, K_rel_array, 'r',label="Non-Wetting Phase")
plt.plot(Snw, K_rel_array1, 'b',label="Wetting Phase")
plt.legend(loc='upper right', ncol=2)
plt.title('Relative K-Phase Saturation', fontsize=18)
plt.xlabel(' Wetting Phase Saturation', fontsize=20)
plt.ylabel('Relative K', fontsize=20)
plt.savefig('plot-pc.png')
plt.show()
"""
"""
book = xlwt.Workbook(encoding="utf-8")

sheet1 = book.add_sheet(extracted_network)

sheet1.write(0, 0, "Sw")

sheet1.write(0, 1, "Krw")
sheet1.write(0, 2, "Krn")
sheet1.write(0, 3, "Sw")
sheet1.write(0, 4, "Pc")
sheet1.write(0, 5, "Snw")
print(K_rel_array)

for i in range(len(Pc)):
    i = i+1
    sheet1.write(i, 4, Pc[i-1])
    sheet1.write(i, 3, Sw[i-1])

    sheet1.write(i, 0, Sw[i-1])
    sheet1.write(i,5,Snw[i-1])

j=0
for n in K_rel_array:
    j = j+1
    sheet1.write(j, 2, n)
k=0
for n in K_rel_array1:
    k = k+1
    sheet1.write(k, 1, n)


book.save("both_throat_pore.xls")


rb = open_workbook('bentheimer-3.18.xls',formatting_info=True)
rs = rb.sheet_by_index(0)

wb = copy(rb)

sheet1 = wb.add_sheet(extracted_network)
sheet1.write(0, 0, "Sw")
sheet1.write(0, 1, "Krw")
sheet1.write(0, 2, "Krn")
sheet1.write(0, 3, "Sw")
sheet1.write(0, 4, "Pc")
sheet1.write(0, 5, "Snw")
print(K_rel_array)

for i in range(len(Pc)):
    i = i+1
    sheet1.write(i, 4, Pc[i-1])
    sheet1.write(i, 3, Sw[i-1])

    sheet1.write(i, 0, Sw[i-1])
    sheet1.write(i,5,Snw[i-1])
j=0
for n in K_rel_array:
    j = j+1
    sheet1.write(j, 2, n)
k=0
for n in K_rel_array1:
    k = k+1
    sheet1.write(k, 1, n)

wb.save("bentheimer-3.18.xls")

"""
