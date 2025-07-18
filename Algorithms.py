import numpy as np
from numpy.linalg import inv
from scipy.stats import beta
from vector_class import Vector
from matplotlib import pyplot as plt

# 所需算法类(找的)

def LeastSquares(A, b):
    transposed = A.T
    inverted = inv(transposed @ A)
    x = (inverted @ transposed) @ b
    return x

def light_vector(magnitude, alpha, beta):
    
    """
    alpha: angle [degs] of light vector with x-axis in x-z plane
    beta: angle [degs] of light vector with x-axis in x-y plane
    """
    alpha = np.deg2rad(alpha)
    beta = np.deg2rad(beta)
    x_component = magnitude * np.cos(alpha) * np.cos(beta)
    y_component = magnitude * np.cos(alpha) * np.sin(beta)
    z_component = - magnitude * np.sin(alpha)
    light_vector = Vector((x_component, y_component, z_component))
    return light_vector