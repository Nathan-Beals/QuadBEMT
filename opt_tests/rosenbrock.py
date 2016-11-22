"""
Solves Rosenbrock's Unconstrained Problem

    min:    100*(x2-x1^2)**2 + (1-x1)^2
    s.t.:   -10 <= xi <= 10, i = 1,2

    f* = 0, x* = [1, 1]
"""


from pyOpt import Optimization
from pyOpt import PSQP
from pyOpt import SLSQP
from pyOpt import SLSQP
from pyOpt import CONMIN
from pyOpt import COBYLA
from pyOpt import SOLVOPT
from pyOpt import KSOPT
from pyOpt import NSGA2
from pyOpt import SDPEN


def objfun(x):

    f = 100*(x[1]-x[0]**2)**2 + (1-x[0])**2
    g = []
    fail = 0

    return f, g, fail


opt_prob = Optimization('Rosenbrock Unconstrained Problem', objfun)
opt_prob.addVar('x1', 'c', lower=-10.0, upper=10.0, value=0.0)
opt_prob.addVar('x2', 'c', lower=-10.0, upper=10.0, value=0.0)
opt_prob.addObj('f')
print opt_prob

# Instantiate optimizer (PSQP) and solve problem
psqp = PSQP()
psqp.setOption('IPRINT', 0)
psqp(opt_prob, sens_type='FD')
print opt_prob.solution(0)

# Instantiate optimizer (SLSQP) and solve problem
slsqp = SLSQP()
slsqp.setOption('IPRINT', -1)
slsqp(opt_prob, sens_type='FD')
print opt_prob.solution(1)

# Instantiate optimizer (CONMIN) and solve problem
conmin = CONMIN()
conmin.setOption('IPRINT', 0)
conmin(opt_prob, sens_type='CS')
print opt_prob.solution(2)

# Instantiate optimizer (COBYLA) and solve problem
cobyla = COBYLA()
cobyla.setOption('IPRINT', 0)
cobyla(opt_prob)
print opt_prob.solution(3)

# Instantiate optimizer (SOLVOPT) and solve problem
solvopt = SOLVOPT()
solvopt.setOption('iprint', -1)
solvopt(opt_prob, sens_type='FD')
print opt_prob.solution(4)

# Instantiate optimizer (KSOPT) and solve problem
ksopt = KSOPT()
ksopt.setOption('IPRINT', 0)
ksopt(opt_prob, sens_type='FD')
print opt_prob.solution(5)

# NSGA2
nsga2 = NSGA2()
nsga2.setOption('PrintOut', 0)
nsga2(opt_prob)
print opt_prob.solution(6)

# SDPEN
sdpen = SDPEN()
sdpen.setOption('iprint', -1)
sdpen(opt_prob)
print opt_prob.solution(7)