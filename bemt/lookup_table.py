import os
import numpy as np
from unit_conversion import deg2rad


def create_table(airfoil_name):
    airfoil_data_dir = 'airfoil_data/' + airfoil_name
    files = os.listdir(airfoil_data_dir)
    reynolds_numbers = []
    alphas = []
    CLs = []
    CDs = []
    i = 0
    for f in files:
        file_path = airfoil_data_dir + '/' + f
        Re = 0
        # Get Reynolds number
        with open(file_path) as fp:
            for i, line in enumerate(fp):
                if i == 8:
                    line = line.replace(" ", "")
                    start_index = line.index('Re=') + 3
                    end_index = line.index('Ncrit') - 1
                    Re = float(line[start_index:end_index+1])
        alpha, CL, CD = np.loadtxt(file_path, dtype=float, skiprows=12, usecols=(0, 1, 2), unpack=True)
        reynolds_numbers += [Re]*len(alpha)
        alphas += list(alpha)
        CLs += list(CL)
        CDs += list(CD)
        i += 1
    reynolds_numbers = np.array(reynolds_numbers)
    alphas = np.array(alphas)
    CLs = np.array(CLs)
    CDs = np.array(CDs)
    return alphas, reynolds_numbers, CLs, CDs