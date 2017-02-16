import os
import numpy as np
from unit_conversion import deg2rad
from scipy.interpolate import griddata
import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
import itertools


def create_table(airfoil_name):
    airfoil_data_dir = 'airfoil_data/' + airfoil_name
    files = os.listdir(airfoil_data_dir)
    reynolds_numbers = np.array([])
    alphas = np.array([])
    CLs = np.array([])
    CDs = np.array([])
    step = 0.1
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

        order = alpha.argsort()
        alpha = alpha[order]
        CL = CL[order]
        CD = CD[order]

        alpha_full = np.arange(-10., 20.1, step)  # Change the step size to reflect the step size polar files
        alpha_pad_low = np.arange(alpha_full[0], alpha[0], step)
        alpha_pad_high = np.arange(alpha[-1]+step, alpha_full[-1]+step, step)
        alpha_table = np.concatenate((alpha_pad_low, alpha, alpha_pad_high))
        coeff_pad_low = np.empty(len(alpha_pad_low))*np.nan
        coeff_pad_high = np.empty(len(alpha_pad_high))*np.nan
        CL_table = np.concatenate((coeff_pad_low, CL, coeff_pad_high))
        CD_table = np.concatenate((coeff_pad_low, CD, coeff_pad_high))

        reynolds_numbers = np.concatenate((reynolds_numbers, np.ones(len(alpha_table))*Re))
        alphas = np.concatenate((alphas, alpha_table))
        CLs = np.concatenate((CLs, CL_table))
        CDs = np.concatenate((CDs, CD_table))
        i += 1
    return alphas, reynolds_numbers, CLs, CDs


def interp_weights(xy, uv, d=2):
    tri = qhull.Delaunay(xy)
    simplex = tri.find_simplex(uv)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uv - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))


# def interpolate(values, vtx, wts):
#     return np.einsum('nj,nj->n', np.take(values, vtx), wts)


def interpolate(t, alpha, Re):
    return griddata(t[0], t[1], (alpha, Re))