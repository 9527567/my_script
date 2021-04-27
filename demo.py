#!/usr/bin/python3
from rdkit import Chem
from rdkit.Chem import AllChem
import os
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit import DataStructs
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

# 读取小分子库


def readNsdf(filename):
    mols_suppl = Chem.SDMolSupplier(filename)
    mols_free = [x for x in mols_suppl if x is not None]
    print(mols_free[0])
    return mols_free


# 读取已知活性小分子化合物，建立参照药效团指纹
def read1sdf(filename):
    mol = Chem.SDMolSupplier(filename)
    if len(mol) > 1:
        raise Exception("can't read more than one molecule for reference!")
    else:
        return mol


# 产生3D构象,如何使用？
def get3D(mol):
    if len(mol) == 1:
        m2 = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        m3 = Chem.RemoveHs(m2)
        return m3
    else:
        result_3D = []
        for x in mol:
            temp = Chem.AddHs(x)
            AllChem.EmbedMolecule(x)
            temp3 = Chem.RemoveHs(temp)
            result.append(temp3)

        return result_3D


# 聚类算法,使用外部库


# 写操作，输出
def writesdf(writename):
    pass


# 获取指纹
def getFinger(mols_suppl, method):
    fps = []
    if method == "top":
        for x in mols_suppl:
            try:
                fps.append(Chem.RDKFingerprint(x))
            except:
                continue
    elif method == "MACCS":
        for x in mols_suppl:
            try:
                fps.append(MACCSkeys.GenMACCSKeys(x))
            except:
                continue
    elif method == "Atoms_pairs":
        for x in mols_suppl:
            try:
                fps.append(Chem.AtomPairs.Pairs.GetAtomPairFingerprint(x))
            except:
                continue
    elif method == "top_torsions":
        for x in mols_suppl:
            try:
                fps.append(Chem.AtomPairs.Torsions.GetTopologicalTorsionFingerprintAsIntVect(x))
            except:
                continue
    else:
        raise Exception("Method is bad!")
    
    return fps


# 计算相似性,相似性方法与指纹方法并不独立，存在依赖关系
def FingerSimilarity(fps, ref_fps, method):
    if method == "Tanimoto_top":
        similarity = [DataStructs.FingerprintSimilarity(ref_fps, i) for i in fps]
    elif method == "Dice_MACCS":
        similarity = [DataStructs.FingerprintSimilarity(ref_fps, i,
                                                        metric=DataStructs.DiceSimilarity) for i in fps]
    elif method == "other":
        raise Exception("support bad!")

    return similarity

#main,独立更多的方法出来

def __main__():
    if len(os.sys.argv) == 1:
        raise Exception("please input molecule library!")
    elif len(os.sys.argv) == 2:
        if os.sys.argv[1] == "--help" or os.sys.argv[1] == "-h":
            print("Using it is easy. chack source code then you can run it!")
        else:
            raise Exception("please input reference molecule!")
    else:
        mols_suppl = readNsdf(os.sys.argv[1])
        print("get", len(mols_suppl), "mol!")
        # fps = getFinger(mols_suppl,method = "MACCS")
        # FingerSimilarity(fps, fps[0],method = "Dice_MACCS")


__main__()
