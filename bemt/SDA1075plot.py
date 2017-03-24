from lookup_table import create_table
from matplotlib import pyplot as plt
import numpy as np


alphas, reynolds_numbers, CLs, CDs, Clmax = create_table('SDA1075_494p')

Re20k_CL = []
Re20k_CD = []
Re20k_alpha = []

Re30k_CL = []
Re30k_CD = []
Re30k_alpha = []

Re40k_CL = []
Re40k_CD = []
Re40k_alpha = []

Re50k_CL = []
Re50k_CD = []
Re50k_alpha = []

Re60k_CL = []
Re60k_CD = []
Re60k_alpha = []

Re100k_CL = []
Re100k_CD = []
Re100k_alpha = []

Re250k_CL = []
Re250k_CD = []
Re250k_alpha = []

for i, re_num in enumerate(reynolds_numbers):
    if alphas[i] >= -5.0:
        if re_num == 20000:
            Re20k_CL.append(CLs[i])
            Re20k_CD.append(CDs[i])
            Re20k_alpha.append(alphas[i])
        elif re_num == 30000:
            Re30k_CL.append(CLs[i])
            Re30k_CD.append(CDs[i])
            Re30k_alpha.append(alphas[i])
        elif re_num == 40000:
            Re40k_CL.append(CLs[i])
            Re40k_CD.append(CDs[i])
            Re40k_alpha.append(alphas[i])
        elif re_num == 50000:
            Re50k_CL.append(CLs[i])
            Re50k_CD.append(CDs[i])
            Re50k_alpha.append(alphas[i])
        elif re_num == 60000:
            Re60k_CL.append(CLs[i])
            Re60k_CD.append(CDs[i])
            Re60k_alpha.append(alphas[i])
        elif re_num == 100000:
            Re100k_CL.append(CLs[i])
            Re100k_CD.append(CDs[i])
            Re100k_alpha.append(alphas[i])
        elif re_num == 250000:
            Re250k_CL.append(CLs[i])
            Re250k_CD.append(CDs[i])
            Re250k_alpha.append(alphas[i])

CD_theory = [0.02 - 0.0216*(alpha*2*np.pi/360) + 0.600*(alpha*2*np.pi/360)**2 for alpha in Re100k_alpha]

# plt.figure(1)
# plt.plot(Re20k_CD, Re20k_CL, Re30k_CD, Re30k_CL, Re40k_CD, Re40k_CL, Re50k_CD, Re50k_CL, Re60k_CD, Re60k_CL, Re100k_CD, Re100k_CL, Re250k_CD, Re250k_CL)
# plt.xlim([0.00, 0.05])
# plt.ylim([0.00, 1.5])
# plt.xlabel("C_D")
# plt.ylabel("C_L")
# plt.legend(["Re=20k", "Re=30k", "Re=40k", "Re=50k", "Re=60k", "Re=100k"], loc='upper left')

plt.figure(2)
plt.plot(Re20k_alpha, Re20k_CL, 'k*-', markerfacecolor='white', markevery=6)
plt.plot(Re40k_alpha, Re40k_CL, 'kv-', markerfacecolor='white', markevery=6)
plt.plot(Re60k_alpha, Re60k_CL, 'ks-', markerfacecolor='white', markevery=6)
plt.plot(Re100k_alpha, Re100k_CL, 'ko-', markerfacecolor='white', markevery=6)
plt.plot(Re250k_alpha, Re250k_CL, 'kD-',markerfacecolor='white', markevery=6)
plt.xlim([-5., 25.])
plt.ylim([-1.0, 2.0])
plt.xlabel(r'$\alpha,\,\mathrm{deg}$', fontsize=18)
plt.ylabel(r'$\mathrm{C}_\mathrm{l}$', fontsize=18)
plt.legend(["Re=20k", "Re=40k", "Re=60k", "Re=100k", "Re=250k"], loc='upper left')
plt.tick_params(axis='both', which='major', labelsize=14)
plt.tick_params(axis='both', which='minor', labelsize=14)
plt.grid()

plt.figure(3)
plt.plot(Re20k_alpha, Re20k_CD, 'k*-', markerfacecolor='white', markevery=6)
plt.plot(Re40k_alpha, Re40k_CD, 'kv-', markerfacecolor='white', markevery=6)
plt.plot(Re60k_alpha, Re60k_CD, 'ks-', markerfacecolor='white', markevery=6)
plt.plot(Re100k_alpha, Re100k_CD, 'ko-', markerfacecolor='white', markevery=6)
plt.plot(Re250k_alpha, Re250k_CD, 'kD-',markerfacecolor='white', markevery=6)
plt.xlim([-5., 25.])
plt.xlabel(r'$\alpha,\,\mathrm{deg}$', fontsize=18)
plt.ylabel(r'$\mathrm{C}_\mathrm{d}$', fontsize=18)
plt.legend(["Re=20k", "Re=40k", "Re=60k", "Re=100k", "Re=250k"], loc='upper left')
plt.tick_params(axis='both', which='major', labelsize=14)
plt.tick_params(axis='both', which='minor', labelsize=14)
plt.grid()

plt.show()

