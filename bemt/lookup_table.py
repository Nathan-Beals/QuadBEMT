import os
import numpy as np
from unit_conversion import deg2rad
from scipy.interpolate import griddata
import math


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
        alpha_full = list(np.arange(-10., 20.1, 0.1))  # Change the step size to reflect the step size polar files
        CL_full = [float('nan')] * len(alpha_full)
        CD_full = [float('nan')] * len(alpha_full)
        alpha, CL, CD = np.loadtxt(file_path, dtype=float, skiprows=12, usecols=(0, 1, 2), unpack=True)
        for i, a1 in enumerate(alpha_full):
            for j, a2 in enumerate(alpha):
                if isclose(a1, a2, abs_tol=0.001):
                    CL_full[i] = CL[j]
                    CD_full[i] = CD[j]
                    break
        reynolds_numbers += [Re]*len(alpha_full)
        alphas += alpha_full
        CLs += CL_full
        CDs += CD_full
        i += 1
    reynolds_numbers = np.array(reynolds_numbers)
    alphas = np.array(alphas)
    CLs = np.array(CLs)
    CDs = np.array(CDs)
    return alphas, reynolds_numbers, CLs, CDs


def interpolate(t, alpha, Re):
    return griddata(t[0], t[1], (alpha, Re))


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)