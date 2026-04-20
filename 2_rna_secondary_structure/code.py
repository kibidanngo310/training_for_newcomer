from typing import List, Tuple, Union
import numpy.typing as npt
import numpy as np

def enumerate_pairs(fastafile: str) -> List[Tuple[int, int]]:
    # 課題 2-1
    seq_list = []

    with open(fastafile, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                seq_list.append(line.strip().upper())

    seq = "".join(seq_list)
    seq = seq.replace("T", "U")

    results = []
    for i in range(0, len(seq)):
        for l in range(i, len(seq)):
            if seq[i] == "A" and seq[l] == "U":
                results.append((i+1, l+1))
            if seq[i] == "U" and seq[l] == "A":
                results.append((i+1, l+1))
            if seq[i] == "G" and seq[l] == "C":
                results.append((i+1, l+1))
            if seq[i] == "C" and seq[l] == "G":
                results.append((i+1, l+1))
    return results

def enumerate_possible_pairs(fastafile: str, min_distance: int=4) -> List[Tuple[int, int]]:
    # 課題 2-2
    seq_list = []

    with open(fastafile, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                seq_list.append(line.strip().upper())

    seq = "".join(seq_list)
    seq = seq.replace("T", "U")

    results = []
    for i in range(0, len(seq)):
        for l in range(i, len(seq)):
            if seq[i] == "A" and seq[l] == "U":
                if l - i >= min_distance:
                    results.append((i+1, l+1))
                
            if seq[i] == "U" and seq[l] == "A":
                if l - i >= min_distance:
                    results.append((i+1, l+1))

            if seq[i] == "G" and seq[l] == "C":
                if l - i >= min_distance:
                    results.append((i+1, l+1))

            if seq[i] == "C" and seq[l] == "G":
                if l - i >= min_distance:
                    results.append((i+1, l+1))
    return results

def enumerate_continuous_pairs(fastafile: str, min_distance: int=4, min_length: int=2) -> List[Tuple[int, int, int]]:
    # 課題 2-3

    return []

def create_dotbracket_notation(fastafile: str, min_distance: int=4, min_length: int=2) -> str:
    # 課題 2-4
    return ""

if __name__ == "__main__":
    filepath = "data/AUCGCCAU.fasta"
    # 課題 2-1
    print(enumerate_pairs(filepath))
    # 課題 2-2
    print(enumerate_possible_pairs(filepath))
    # 課題 2-3
    print(enumerate_continuous_pairs(filepath, 2))
    # 課題 2-4
    print(create_dotbracket_notation(filepath, 2))


