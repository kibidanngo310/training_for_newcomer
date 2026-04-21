from typing import List, Union
import numpy.typing as npt
import numpy as np

def calc_chain_center(pdbfile: str, chain: str) -> Union[npt.NDArray[np.float_], List[float]]:
    # 課題 3-2
    return np.zeros(3)

def calc_residue_center(pdbfile: str, chain: str, resname: str) -> Union[npt.NDArray[np.float_], List[float]]:
    # 課題 3-3
    return np.zeros(3)

def calc_min_distance(pdbfile: str, chain: str, resname: str) -> float:
    # 課題 3-4
    #0.406 (827 to 827 atoms)
    return 0.0



if __name__ == "__main__":
    filepath = "data/1buw.pdb"
    chain = "A"
    resname = "HEM"
    # 課題 3-2
    print(calc_chain_center(filepath, chain))
    # 課題 3-3
    print(calc_residue_center(filepath, chain, resname))
    # 課題 3-4
    print(calc_min_distance(filepath, chain, resname))
