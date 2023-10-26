import numpy as np
import scipy


def euclid_to_euler(x, y, z):
    return np.arctan2((x**2 + y**2) ** 0.5, z), np.arctan2(y, x)


def qlm(l, m, x, y, z):
    """
    球面ベクトルの球面調和関数展開
    """
    theta, phi = euclid_to_euler(x, y, z)
    # print(theta,phi)
    return scipy.special.sph_harm(m, l, phi, theta, out=None).mean()


def boo(l, x, y, z):
    """
    結合配向秩序指標
    """
    s = 0
    for m in range(-l, l + 1):
        q = qlm(l, m, x, y, z)
        s += q * q.conj()
    return (s * 4 * np.pi / (2 * l + 1)) ** 0.5
