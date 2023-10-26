import numpy as np
import scipy
import matplotlib.pyplot as plt


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


if __name__ == "__main__":
    # 単純立方格子。隣接ベクトルは6本。
    scl = np.array(
        [
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    # 体心立方格子
    bcc = (
        np.array(
            [
                [-1, -1, -1],
                [-1, -1, +1],
                [-1, +1, -1],
                [-1, +1, +1],
                [+1, -1, -1],
                [+1, -1, +1],
                [+1, +1, -1],
                [+1, +1, +1],
                # Steinhardtとあわせるために、第二隣接までとる。
                [-2.0, 0.0, 0.0],
                [0.0, -2.0, 0.0],
                [0.0, 0.0, -2.0],
                [2.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, 0.0, 2.0],
            ]
        )
        / 3**0.5
    )
    # ダイヤモンド格子
    dia = np.array([[-1, -1, -1], [+1, -1, -1], [+1, -1, +1], [+1, +1, -1]]) / 3**0.5
    # 面心立方格子
    fcc = (
        np.array(
            [
                [-1, -1, 0],
                [-1, +1, 0],
                [+1, -1, 0],
                [+1, +1, 0],
                [-1, 0, -1],
                [-1, 0, +1],
                [+1, 0, -1],
                [+1, 0, +1],
                [0, -1, -1],
                [0, -1, +1],
                [0, +1, -1],
                [0, +1, +1],
            ]
        )
        / 2**0.5
    )

    samples = {"BCC": bcc, "SCL": scl, "Diamond": dia, "FCC": fcc}
    for label, vecs in samples.items():
        x, y, z = vecs.T

        Y = []
        for l in range(13):
            Y.append(boo(l, x, y, z))

        Y = np.array(Y)

        plt.plot(Y.real, label=label)

    plt.legend()
    plt.show()
