import numpy as np
from numpy.linalg import inv


def LeastSquares(A, b):
    # 简单的最小二乘法
    # 假定A列满秩,则使用法方程法以求解。
    transposed = A.T
    inverted = inv(transposed @ A)
    x = (inverted @ transposed) @ b
    return x

def lsq(A,b,rcond):
    # 基于奇异值分解的最小二乘法,带有截断:
    # 小于最大奇异值rcond倍则不予考虑，设为0.
    [U,S,V] = np.linalg.svd(A)
    n,_ = U.shape
    m,_ = V.shape
    Sver = np.zeros((m,n))
    for i in range(0,min(n,m)):
        if S[i] < rcond * S[0]:
            Sver[i,i] = 0
        else:
            Sver[i,i] = 1/S[i]
    return V.T @ Sver @ U.T @ b