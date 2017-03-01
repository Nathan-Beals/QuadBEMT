from pyOpt import Optimization
from pyOpt import CONMIN


def obj_fun(x):

    f = x[0]**2 - 5*x[0] + x[1]**2 - 5*x[1] + 2*x[2]**2 - 21*x[2] + x[3]**2 + 7*x[3] + 50

    g = [0] * 3
    g[0] = x[0]**2 + x[0] + x[1]**2 - x[1] + x[2]**2 + x[2] + x[3]**2 - x[3] - 8
    g[1] = x[1]**2 - x[0] + 2*x[1]**2 + x[2]**2 + 2*x[3]**2 - x[3] - 10
    g[2] = 2*x[0]**2 + 2*x[0] + x[1]**2 - x[1] + x[2]**2 - x[3] - 5

    fail = 0

    return f, g, fail


xinit = [1.0, 1.0, 1.0, 1.0]


opt_prob = Optimization('Rosen-Suzuki Constrained Minimization', obj_fun)
opt_prob.addVarGroup('x', 4, 'c', value=xinit)
opt_prob.addConGroup('g', 3, 'i')
opt_prob.addObj('f')
print opt_prob

conmin = CONMIN()
conmin.setOption('IPRINT', 1)
conmin(opt_prob)
print opt_prob.solution(0)
