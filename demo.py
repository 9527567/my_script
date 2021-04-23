from rdkit import Chem
from rdkit.Chem import AllChem
import os
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs


# 读取小分子库
def readNsdf(filename):
    mols_suppl = Chem.SDMolSupplier(filename)
    return mols_suppl


# 读取已知活性小分子化合物，建立参照药效团指纹
def read1sdf(filename):
    mol = Chem.SDMolSupplier(filename)
    if len(mol) > 1:
        raise Exception("can't read more than one molecule for reference!")
    else:
        return mol


def writesdf(writename):
    pass


# 获取指纹
def getFinger(mols_suppl):
    fps = []
    for x in mols_suppl:
        try:
            fps.append(MACCSkeys.GenMACCSKeys(x))
        except:
            continue
    return fps


# 计算相似性,多种指纹计算方式
def FingerSimilarity(fps):
    temp = [DataStructs.FingerprintSimilarity(fps[0], i,
                                              metric=DataStructs.DiceSimilarity)for i in fps]
    return temp


def main():
    if len(os.sys.argv) == 1:
        raise Exception("please input molecule library!")
    elif len(os.sys.argv) == 2:
        if os.sys.argv[2] == "--help" or os.sys.argv[2] == "-h":
            print("Using it is easy.chack source code then you can run it!")
        else:
            raise Exception("please input reference molecule!")
    else:
        mols_suppl = readNsdf(os.sys.argv[1])
        print("get", len(mols_suppl), "mol!")
        fps = getFinger(mols_suppl)
        FingerSimilarity(fps)


main()
