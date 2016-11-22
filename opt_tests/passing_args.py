from pyOpt import Optimization
from pyOpt import SLSQP


def objfun(x, **kwargs):

    a1 = kwargs['a12'][0]
    a2 = kwargs['a12'][1]
    a3 = kwargs['a3']

    f = a1*(x[1]-x[0]**2.)**2. + (a2-x[0])**2.
    g = [0.0]*2
    g[0] = x[0]**2. + x[1]**2.0 - a3

    fail = 0
    return f, g, fail

# Instantiate optimization problem
opt_prob = Optimization('Rosenbrock Constrained Problem', objfun)
opt_prob.addVar('x1', 'c', lower=0.0, upper=1.0, value=0.5)
opt_prob.addVar('x2', 'c', lower=0.0, upper=1.0, value=0.5)
opt_prob.addObj('f')
opt_prob.addCon('g1', 'i')
print opt_prob

# Arguments to pass to objfun
a1 = 100.0
a2 = 1.0
a3 = 1.0

# Instantiate optimizer (SLSQP) and solve problem
slsqp = SLSQP()
slsqp.setOption('IPRINT', -1)
slsqp(opt_prob, sens_type='FD', a12=[a1, a2], a3=a3)
print opt_prob.solution(0)