import numpy as np


def rpm2radps(n):
    return n * 2 * np.pi / 60


def radps2rpm(n):
    return n * 60 / 2 / np.pi


def rad2deg(n):
    return n * 360 / 2 / np.pi


def deg2rad(n):
    return n * 2 * np.pi / 360


def in2m(n):
    return n * 0.0254

