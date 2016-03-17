import scipy as sp
from scipy import interpolate


x = sp.array([1, 2, 3])
y = sp.array([70000, 80000, 90000])
z = sp.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

f = interpolate.interp2d(x, y, z)

print f(1.5, 75000)
