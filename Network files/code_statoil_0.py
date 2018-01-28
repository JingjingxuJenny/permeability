import numpy as np
import scipy as sp
import OpenPNM as op
import sys
import my_models
print(sys.path)


extracted_network = 'mtsimon-2.55'

pn = op.Utilities.IO.Statoil.load(path=r'./', prefix=extracted_network)
pn.name = extracted_network
op.export_data(network=pn, filename=extracted_network)
print(op,"before")

print(pn.check_network_health()['trim_pores'],"check health")
pn.trim(pn.check_network_health()['trim_pores'])



print(pn,"after")
geom = op.Geometry.GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)
water = op.Phases.Water(network=pn)
phys = op.Physics.GenericPhysics(network=pn, phase=water, geometry=geom)
for item in pn.props():
	if item not in ['throat.conns', 'pore.coords']:
		geom.update({item: pn.pop(item)})
geom['pore.diameter'] = 2*geom['pore.radius']
geom.models.add(propname='pore.area',model=op.Geometry.models.pore_area.spherical)
geom['throat.diameter'] = 2*geom['throat.radius']
geom.models.add(propname='throat.area',model=op.Geometry.models.throat_area.cylinder)
phys.models.add(propname='throat.hydraulic_conductance',model=my_models.noncircular_capillary)
#phys.models.add(propname='throat.hydraulic_conductance',model=op.Physics.models.hydraulic_conductance.hagen_poiseuille)


flow = op.Algorithms.StokesFlow(network=pn, phase=water)
flow.set_boundary_conditions(pores=pn.pores('inlets'), bctype='Dirichlet', bcvalue=200000)
flow.set_boundary_conditions(pores=pn.pores('outlets'), bctype='Dirichlet', bcvalue=100000)
flow.run()
flow.return_results()

K = flow.calc_eff_permeability()
mypermfactor = (1.013249966e+12)*1000
myK = K*mypermfactor
print('Abs. Permeability (mD): ', myK)




"""
## create a new label array for triangle

#print(np.size(geom["throat.shape_factor"] ))
Ts1 = geom["throat.shape_factor"]<0.0481
#print(pn.labels())
pn['throat.triangle'] = Ts1
#print(pn.labels())
#print(np.size(pn['throat.triangle']))
print("throat satisfy triangle ",sp.where(pn['throat.triangle'])[0],"number :",np.size(sp.where(pn['throat.triangle'])[0]))

## create a new label array for rectangle
Ts2 = geom["throat.shape_factor"]>0.0481
Ts3 = geom["throat.shape_factor"]<0.0625
Ts4 = Ts2 & Ts3
pn['throat.rectangle'] = Ts4
print("throat satisfy rectangle ",sp.where(pn['throat.rectangle'])[0],"number:",np.size(sp.where(pn['throat.rectangle'])[0]))


## create a new label array for circle
Ts5 = geom["throat.shape_factor"]>0.0625
pn['throat.circle'] = Ts5
print("throat satisfy circle ",sp.where(pn['throat.circle'])[0],"number:",np.size(sp.where(pn['throat.circle'])[0]))

## since almost the  throat shape factor located in the range of triangle,
## the factor 0.6 is used for all shape factor fo the sake of convenience

print(pn['throat.shape_factor'])
print(pn["throat.triangle"])
# create a new label array for triangle



#print(np.size(geom["throat.shape_factor"] ))
Ts1 = geom["pore.shape_factor"]<0.0481
#print(pn.labels())
pn['pore.triangle'] = Ts1
#print(pn.labels())
#print(np.size(pn['throat.triangle']))
print("pore satisfy triangle ",sp.where(pn['pore.triangle'])[0],"number :",np.size(sp.where(pn['pore.triangle'])[0]))

## create a new label array for rectangle
Ts2 = geom["pore.shape_factor"]>0.0481
Ts3 = geom["pore.shape_factor"]<0.0625
Ts4 = Ts2 & Ts3
pn['pore.rectangle'] = Ts4
print("pore satisfy rectangle ",sp.where(pn['pore.rectangle'])[0],"number:",np.size(sp.where(pn['pore.rectangle'])[0]))


## create a new label array for circle
Ts5 = geom["pore.shape_factor"]>0.0625
pn['pore.circle'] = Ts5
print("pore satisfy circle ",sp.where(pn['pore.circle'])[0],"number:",np.size(sp.where(pn['pore.circle'])[0]))

## since almost the  throat shape factor located in the range of triangle,
## the factor 0.6 is used for all shape factor fo the sake of convenience

print(pn['pore.shape_factor'])
print(pn["pore.triangle"])
"""
