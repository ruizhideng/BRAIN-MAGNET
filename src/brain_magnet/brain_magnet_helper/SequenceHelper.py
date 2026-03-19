"""
Sequence helper functions adapted from DeepSTARR.

Source:
- DeepSTARR `helper/SequenceHelper.py`
  https://github.com/bernardo-de-almeida/DeepSTARR/tree/main/DeepSTARR/Neural_Network_DNA_Demo/helper
"""


def parse_alpha_to_seq(sequence: str) -> "np.ndarray":
    import numpy as np
    output = np.arange(len(sequence))
    for i in range(0, len(sequence)):
        snippet = sequence[i]
        if snippet == "A":
            output[i] = 0
        elif snippet == "C":
            output[i] = 1
        elif snippet == "G":
            output[i] = 2
        elif snippet == "T":
            output[i] = 3
        elif snippet == "N":
            output[i] = -1
        else:
            raise AssertionError("Cannot handle snippet: " + snippet)
    return output


def to_categorical(y: "np.ndarray", nb_classes: int = None) -> "np.ndarray":
    """Convert class vector to binary class matrix."""
    import numpy as np
    y = np.asarray(y, dtype="int32")
    if nb_classes is None:
        nb_classes = np.max(y) + 1
    Y = np.zeros((len(y), nb_classes))
    for i in range(len(y)):
        if y[i] != -1:
            Y[i, y[i]] = 1.0
    return Y


def do_one_hot_encoding(sequence: "np.ndarray", seq_length: int, f=parse_alpha_to_seq) -> "np.ndarray":
    import numpy as np
    X = np.zeros((sequence.shape[0], seq_length, 4))
    for idx in range(0, len(sequence)):
        X[idx] = to_categorical(f(sequence[idx]), 4)
    return X


def generate_complementary_sequence(sequence: str) -> str:
    comp_seq = []
    for b in sequence:
        if b == "A":
            comp_seq.append("T")
        elif b == "T":
            comp_seq.append("A")
        elif b == "C":
            comp_seq.append("G")
        elif b == "G":
            comp_seq.append("C")
        elif b == "N":
            comp_seq.append("N")
        else:
            raise ValueError("Cannot convert base {0} to complement base!".format(b))
    return "".join(comp_seq)
