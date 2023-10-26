#!/usr/bin/env python

import sys

import numpy as np


def read_gro(file):
    """
    gromacsの.groファイルを読みこむ。

    あとで出力する場合にそなえ、できるだけデータをそのままの形で保持する。
    """

    # 無限ループ
    while True:
        frame = {
            "resi_id": [],
            "residue": [],
            "atom": [],
            "atom_id": [],
            "position": [],
        }

        title = file.readline()
        # 終了判定。1文字も読めない時はファイルの終わり。
        if len(title) == 0:
            return
        n_atom = int(file.readline())
        for i in range(n_atom):
            line = file.readline()
            residue_id = int(line[0:5])
            residue = line[5:10].strip()
            atom = line[10:15].strip()
            atom_id = int(line[15:20])
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            # 速度は省略

            frame["resi_id"].append(residue_id)
            frame["residue"].append(residue)
            frame["atom"].append(atom)
            frame["atom_id"].append(atom_id)
            frame["position"].append([x, y, z])

        cell = [float(x) for x in file.readline().split()]

        # numpy形式に変換しておく。
        frame["resi_id"] = np.array(frame["resi_id"])
        frame["residue"] = np.array(frame["residue"])
        frame["atom"] = np.array(frame["atom"])
        frame["atom_id"] = np.array(frame["atom_id"])
        frame["position"] = np.array(frame["position"])

        # cellは行列の形にしておく。
        if len(cell) == 3:
            # 直方体セルの場合
            cell = np.diag(cell)
        else:
            # 9パラメータで指定される場合は、順番がややこしい。
            # v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
            x = [cell[0], cell[5], cell[7]]
            y = [cell[3], cell[1], cell[8]]
            z = [cell[4], cell[6], cell[2]]
            cell = np.array([x, y, z])

        frame["cell"] = cell
        # returnの代わりにyieldを使うと、繰り返し(iterator)にできる。
        yield frame


if __name__ == "__main__":
    for frame in read_gro(sys.stdin):
        print(frame)
