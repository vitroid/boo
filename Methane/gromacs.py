#!/usr/bin/env python

import sys
from collections import defaultdict

import numpy as np


def read_gro(file):
    """
    `read_gro` 関数は GROMACS .gro ファイルを読み取り、軌跡の各フレームを辞書として生成します。

    Args:
      file: `file` パラメータは、読み込む GRO ファイルを表すファイル オブジェクトです。 `read_gro` 関数に渡す前に、読み取りモードで開く必要があります。

    Returns:
      関数 `read_gro` は、ジェネレーターを作成するために `yield`
    ステートメントを使用しています。この関数は、値を返す代わりに、ループの反復ごとに「frame」辞書を生成します。これにより、呼び出し元はすべてのフレームを一度に返すのではなく、フレームを 1
    つずつ繰り返すことができます。
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
        # frame["resi_id"] = np.array(frame["resi_id"])
        # frame["residue"] = np.array(frame["residue"])
        # frame["atom"] = np.array(frame["atom"])
        # frame["atom_id"] = np.array(frame["atom_id"])
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


def write_gro(frame, file, remark="Written by write_gro"):
    """
    `write_gro` 関数は、分子動力学シミュレーションで一般的に使用される GRO 形式でフレームをファイルに書き込みます。

    Args:
      frame: 「frame」パラメータは、次のキーを含む辞書です。
      file: `file` パラメータは、フレームが書き込まれるファイル オブジェクトです。 `write_gro` 関数を呼び出す前に、書き込みモードで開く必要があります。
      remark: 「remark」パラメータは、出力ファイルの最初の行として書き込まれるメッセージまたはコメントを表す文字列です。これはオプションであり、デフォルト値は「Written by
    write_gro」です。. Defaults to Written by write_gro
    """

    # frameを解体する
    residue_id = frame["resi_id"]
    residue = frame["residue"]
    atom = frame["atom"]
    atom_id = frame["atom_id"]
    positions = frame["position"]

    # 1行目はメッセージ行
    print(remark, file=file)
    # 2行目は原子数
    Natom = len(positions)
    print(Natom, file=file)
    # 原子もそのまま
    for i in range(Natom):
        ri = residue_id[i]
        r = residue[i]
        a = atom[i]
        ai = atom_id[i]
        pos = positions[i]
        print(
            f"{ri:5d}{r:5s}{a:>5s}{ai:5d}{pos[0]:8.3f}" f"{pos[1]:8.3f}{pos[2]:8.3f}",
            file=file,
        )
    # セルは、直方体とそれ以外で書き方が違う
    cell = frame["cell"]
    if cell[1, 0] == 0:
        print(cell[0, 0], cell[1, 1], cell[2, 2], file=file)
    else:
        print(
            cell[0, 0],
            cell[1, 1],
            cell[2, 2],
            cell[1, 0],
            cell[2, 0],
            cell[0, 1],
            cell[2, 1],
            cell[0, 2],
            cell[1, 2],
            file=file,
        )


def decompose(frame):
    """
    関数「decompose」は原子位置のフレームを取得し、それらを残基 ID に基づいて分子に分離します。

    Args:
      frame: 「frame」パラメータは、分子システム内の原子に関する情報を含む辞書です。次のキーがあります。

    Returns:
      キーが残基の名前であり、値が分子のリストである辞書。各分子は、原子と位置のペアのリストとして表されます。
    """
    Natom = frame["position"].shape[0]
    molecules = defaultdict(list)
    lastres = -1
    molecule = []
    for i in range(Natom):
        resi_id = frame["resi_id"][i]
        if resi_id != lastres:
            if len(molecule) > 0:
                molecules[resname].append(molecule)
            molecule = []
            lastres = resi_id
        resname = frame["residue"][i]
        atom = frame["atom"][i]
        atom_id = frame["atom_id"][i]
        position = frame["position"][i].copy()
        molecule.append([atom, position])
    if len(molecule) > 0:
        molecules[resname].append(molecule)
    return molecules


def compose(mols, cell):
    """
    「compose」関数は、分子とその位置、およびセル サイズの辞書を受け取り、残基 ID、残基名、原子名、原子 ID、位置、およびセル サイズを含む辞書を返します。

    Args:
      mols: `mols` パラメータは、キーが分子の名前、値が原子のリストである辞書です。各原子は、原子名とその位置を含むタプルとして表されます。
      cell: 「セル」パラメータは、分子が構成されるセルの寸法を表します。これは、セルの長さ、幅、高さをそれぞれ表す 3 つの値を含むリストまたはタプルです。

    Returns:
      次のキーと値を含む辞書:
    """
    resi_id = []
    residue = []
    atom = []
    atom_id = []
    position = []
    aid = 0
    rid = 0
    for name, mollist in mols.items():
        for mol in mollist:
            rid += 1
            for a in mol:
                aid += 1
                resi_id.append(rid)
                residue.append(name)
                atom.append(a[0])
                atom_id.append(aid)
                position.append(a[1])
    return {
        "resi_id": resi_id,
        "residue": residue,
        "atom": atom,
        "atom_id": atom_id,
        "position": position,
        "cell": cell,
    }


if __name__ == "__main__":
    for frame in read_gro(sys.stdin):
        print(frame)
