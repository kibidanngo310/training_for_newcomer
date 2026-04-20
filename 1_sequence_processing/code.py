from typing import List, Union
import numpy.typing as npt
import numpy as np

def base_count(fastafile: str) -> List[int]:
    # 課題 1-1

    counts = {'A':0, 'T':0, 'G':0, 'C':0}
    with open(fastafile, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            line = line.strip().upper()

            for base in line:
                if base in counts:
                    counts[base] +=1
    return [counts['A'], counts['T'], counts['G'], counts['C']] # A, T, G, C

def gen_rev_comp_seq(fastafile: str) -> str:
    # 課題 1-2
    comp = {'A':'T','T':'A','G':'C','C':'G'}
    seq_list = []

    with open(fastafile, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq_list.append(line.strip().upper())

    seq = "".join(seq_list)

    reverse = ""
    for base in reversed(seq):
        if base in comp:
            reverse += comp[base]
        else:
            reverse += base
    return reverse

def calc_gc_content(fastafile: str, window: int=1000, step: int=300) -> Union[npt.NDArray[np.float64], List[float]]:
    # 課題 1-3
    # 値を出力するところまで。matplotlibを使う部分は別途実装してください。
    seq_list = []

    with open(fastafile, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                seq_list.append(line.strip().upper())

    seq = "".join(seq_list)
    results = []
    for i in range(0, len(seq) - window + 1, step):
        sub_seq = seq[i: i + window]
        counts = sub_seq.count('G') + sub_seq.count('C')
        results.append((counts/window) * 100)
    return results

def search_motif(fastafile: str, motif: str) -> List[str]:
    # 課題 1-4
    seq_list = []

    with open(fastafile, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                seq_list.append(line.strip().upper())

    seq = "".join(seq_list)
    revseq = gen_rev_comp_seq(fastafile)

    results: List[str] = [] 

    for i in range(0, len(seq) - len(motif) + 1):
        if seq[i :i + len(motif)] == motif:
            results.append(f"F{i + 1}")
    
    for i in range(0, len(revseq) - len(motif) + 1):
        if revseq[i :i + len(motif)] == motif:
            results.append(f"R{len(revseq) - i}")

    return results

def translate(fastafile: str) -> List[str]:
    # 課題 1-5
    codon_table = {
        'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
        'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
        'TAT':'Y', 'TAC':'Y', 'TAA':'_', 'TAG':'_',
        'TGT':'C', 'TGC':'C', 'TGA':'_', 'TGG':'W',

        'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
        'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
        'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
        'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',

        'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
        'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
        'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
        'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',

        'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
        'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
        'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
        'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
    }

    seq_list = []

    with open(fastafile, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                seq_list.append(line.strip().upper())

    seq = "".join(seq_list)
    revseq = gen_rev_comp_seq(fastafile)

    results = []

    for s in [seq, revseq]:
        for frame in range(3):

            protein = ""
            for i in range(frame, len(s) - 2, 3):
                codon = s[i:i+3]
                protein += codon_table.get(codon, "")

            current = ""
            in_protein = False

            for aa in protein:
                if aa == "M":
                    current = "M"
                    in_protein = True
                elif aa == "_":
                    if in_protein:
                        current += "_"
                        results.append(current)
                        in_protein = False
                elif in_protein:
                    current += aa

            if in_protein:
                results.append(current)

    return []

if __name__ == "__main__":
    filepath = "data/NT_113952.1.fasta"
    # 課題 1-1
    print(base_count(filepath))
    # 課題 1-2
    print(gen_rev_comp_seq(filepath))
    # 課題 1-3
    print(calc_gc_content(filepath))
    # 課題 1-4
    print(search_motif(filepath, "ATG"))
    # 課題 1-5
    print(translate(filepath))
