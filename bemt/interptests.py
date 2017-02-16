from lookup_table import create_table, interpolate, interpolate_griddata, interp_weights
from scipy.interpolate import griddata, interp2d
import numpy as np
import timeit


def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped


airfoil = 'SDA1075_494p'
alpha, Re, CL, CD = create_table(airfoil)
xy = (alpha, Re)
u = np.array([5.0, 6.0])
v = np.array([80000, 85000])

print "griddata val = " + str(griddata(xy, CL, (u,v)))
wrapped = wrapper(griddata, (alpha, Re), CL, (5.0, 80000))
print timeit.timeit(wrapped, number=100)

# vtx, wts = interp_weights(xy, uv)
# print "interp val = " + str(interpolate(CL, vtx, wts))


