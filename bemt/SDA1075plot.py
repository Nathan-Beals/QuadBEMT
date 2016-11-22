from lookup_table import create_table
from matplotlib import pyplot as plt


alphas, reynolds_numbers, CLs, CDs = create_table('SDA1075')

Re40k_CL = []
Re40k_CD = []
Re60k_CL = []
Re60k_CD = []
Re100k_CL = []
Re100k_CD = []
for i, re_num in enumerate(reynolds_numbers):
    if re_num == 40000:
        Re40k_CL.append(CLs[i])
        Re40k_CD.append(CDs[i])
    elif re_num == 60000:
        Re60k_CL.append(CLs[i])
        Re60k_CD.append(CDs[i])
    elif re_num == 100000:
        Re100k_CL.append(CLs[i])
        Re100k_CD.append(CDs[i])

plt.plot(Re40k_CD, Re40k_CL, Re60k_CD, Re60k_CL, Re100k_CD, Re100k_CL)
plt.xlim([0.00, 0.05])
plt.ylim([0.00, 1.5])
plt.xlabel("C_D")
plt.ylabel("C_L")
plt.legend(["Re=40000", "Re=60000", "Re=100000"], loc='upper left')
plt.show()