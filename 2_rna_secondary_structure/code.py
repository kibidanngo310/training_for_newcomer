from typing import List, Tuple, Union
import numpy.typing as npt
import numpy as np

def pair(s1: str, s2: str) -> bool:
    result = False
    if s1 == "A" and s2 == "U":
        result = True
    if s1 == "U" and s2 == "A":
        result = True
    if s1 == "G" and s2 == "C":
        result = True
    if s1 == "C" and s2 == "G":
        result = True
    return result

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
            if pair(seq[i], seq[l]):
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
            if pair(seq[i], seq[l]):
                if (l - i >= min_distance):
                    results.append((i+1, l+1))
    return results

def enumerate_continuous_pairs(fastafile: str, min_distance: int=4, min_length: int=2) -> List[Tuple[int, int, int]]:
    # 課題 2-3
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
            count = 0
            sub_i = i
            sub_l = l
            while pair(seq[sub_i], seq[sub_l]) and sub_l - sub_i >= min_distance:
                sub_i += 1
                sub_l -= 1
                count += 1
            if count >= min_length:
                results.append((i+1, l+1, count))

    return results

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


