import scipy as sp
def mason_model(network, phase, physics, f=0.6667, **kwargs):
    Dt = network['throat.diameter']
    theta = phase['throat.contact_angle']
    sigma = phase['throat.surface_tension']
    Pc = 4*sigma*sp.cos(f*sp.deg2rad(theta))/Dt
    return Pc[network.throats(physics.name)]
