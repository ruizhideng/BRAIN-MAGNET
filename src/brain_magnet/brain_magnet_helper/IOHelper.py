"""
IO helper functions adapted from DeepSTARR.

Source:
- DeepSTARR `helper/IOHelper.py`
  https://github.com/bernardo-de-almeida/DeepSTARR/tree/main/DeepSTARR/Neural_Network_DNA_Demo/helper
"""

import gzip
import math
from typing import Dict, Optional, Union


def get_fastas_from_file(
    fasta_path: str,
    as_dict: bool = False,
    uppercase: bool = False,
    stop_at: Optional[int] = None,
):  # type: ignore[no-untyped-def]
    """
    Read FASTA into either a dict or a pandas.DataFrame-like structure.

    Pandas is imported lazily so the core package can be used without pandas.
    """
    fastas = []
    seq = None
    header = None

    opener = gzip.open if fasta_path.endswith(".gz") else open
    with opener(fasta_path, "rt") as f:
        for r in f:
            r = r.strip()
            if not r:
                continue
            if r.startswith(">"):
                if seq is not None and header is not None:
                    fastas.append([header, seq])
                    if stop_at is not None and len(fastas) >= stop_at:
                        break
                seq = ""
                header = r[1:]
            else:
                if seq is None:
                    seq = r.upper() if uppercase else r
                else:
                    seq += r.upper() if uppercase else r

    # append last fasta read by method
    if stop_at is not None:
        if len(fastas) < stop_at and header is not None and seq is not None:
            fastas.append([header, seq])
    else:
        if header is not None and seq is not None:
            fastas.append([header, seq])

    if as_dict:
        return {h: s for h, s in fastas}

    import pandas as pd  # type: ignore

    return pd.DataFrame({"location": [e[0] for e in fastas], "sequence": [e[1] for e in fastas]})


def get_padded_sequences(fasta_file: str):
    fasta = get_fastas_from_file(fasta_file)
    max_length = max(len(x) for x in fasta.sequence)
    padded_sequences = []
    for seq in fasta.sequence:
        diff = max_length - len(seq)
        n_seq = (math.floor(diff / 2) * "N") + seq + (math.ceil(diff / 2) * "N")
        padded_sequences.append(n_seq)
    fasta.sequence = padded_sequences
    return fasta
