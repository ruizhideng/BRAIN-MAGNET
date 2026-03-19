from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Sequence, Tuple, TypeVar, Union


class EnhancerRecord:
    def __init__(self, chrom, start, end, activities):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.activities = activities  # type: List[float]

    @property
    def fasta_header(self):
        return f"{self.chrom}:{self.start}-{self.end}"

    @property
    def reversed_fasta_header(self):
        # Matches the historical format observed in existing training FASTAs:
        # `>107931776-107932306 Reversed:`
        return f"{self.start}-{self.end} Reversed:"


def read_enhancer_activity_tsv(path):  # type: (Union[str, Path]) -> Tuple[List[EnhancerRecord], List[str]]
    path = Path(path)
    records = []  # type: List[EnhancerRecord]
    activity_names = []  # type: List[str]

    with path.open() as f:
        header = f.readline().strip("\n").split("\t")
        if len(header) < 4:
            raise ValueError(
                "Expected a tab-separated file with at least 4 columns: "
                "first 3 are BED (chrom, start, end) followed by >=1 activity columns."
            )

        activity_names = header[3:]

        for line_no, line in enumerate(f, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            try:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                activities = [float(x) for x in parts[3:]]
            except Exception as e:  # noqa: BLE001 (surface parsing context)
                raise ValueError(f"Failed parsing {path}:{line_no}: {line}") from e

            if len(activities) != len(activity_names):
                raise ValueError(
                    "Activity column count mismatch at {p}:{line_no}. "
                    "Header has {h} activity columns, row has {r}.".format(
                        p=path, line_no=line_no, h=len(activity_names), r=len(activities)
                    )
                )

            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}" if chrom else chrom
            if start < 0 or end <= start:
                raise ValueError(f"Invalid interval at {path}:{line_no}: start={start}, end={end}")

            records.append(
                EnhancerRecord(
                    chrom=chrom,
                    start=start,
                    end=end,
                    activities=activities,
                )
            )

    return records, activity_names


_RC_TABLE = str.maketrans(
    {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "N": "N",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a",
        "n": "n",
    }
)


def reverse_complement(seq):  # type: (str) -> str
    return seq.translate(_RC_TABLE)[::-1]


def iter_fasta(
    records: Sequence[EnhancerRecord],
    sequences: Sequence[str],
    *,
    include_reversed: bool,
):  # type: (...) -> Iterator[Tuple[str, str]]
    if len(records) != len(sequences):
        raise ValueError(f"records/sequences length mismatch: {len(records)} vs {len(sequences)}")

    for rec, seq in zip(records, sequences):
        yield rec.fasta_header, seq

    if include_reversed:
        for rec, seq in zip(records, sequences):
            yield rec.reversed_fasta_header, reverse_complement(seq)


def write_fasta(pairs, out_fasta):  # type: (Iterable[Tuple[str, str]], Union[str, Path]) -> None
    out_fasta = Path(out_fasta)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with out_fasta.open("w") as f:
        for header, seq in pairs:
            f.write(f">{header}\n")
            f.write(f"{seq}\n")


def write_activity_tsv(
    records: Sequence[EnhancerRecord],
    out_tsv,
    *,
    include_reversed: bool,
    activity_names,  # type: Sequence[str]
    activity_indices=None,  # type: Optional[Sequence[int]]
) -> None:
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    if activity_indices is None:
        activity_indices = list(range(len(activity_names)))

    with out_tsv.open("w") as f:
        f.write("\t".join([activity_names[i] for i in activity_indices]) + "\n")
        for r in records:
            f.write("\t".join([str(r.activities[i]) for i in activity_indices]) + "\n")
        if include_reversed:
            for r in records:
                f.write("\t".join([str(r.activities[i]) for i in activity_indices]) + "\n")


def fetch_genome_sequences(
    genome_fasta,
    records: Sequence[EnhancerRecord],
):  # type: (...) -> List[str]
    """
    Fetch sequences using 0-based, end-exclusive intervals: [start, end).

    This matches common BED conventions and yields length = end - start.
    Requires ``pyfaidx`` (install with ``pip install pyfaidx`` or via conda/bioconda).
    For maximum portability on HPC clusters, prefer providing ``sequences_fasta``
    (see :func:`load_sequences_from_fasta`) instead.
    """
    genome_fasta = Path(genome_fasta)
    try:
        from pyfaidx import Fasta  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(
            "Genome-based sequence fetching requires `pyfaidx`.\n"
            "Either install it in your environment (if possible), or provide `--sequences-fasta` instead."
        ) from e

    genome = Fasta(str(genome_fasta), as_raw=True, sequence_always_upper=True)
    seqs = []  # type: List[str]
    for r in records:
        try:
            seq = genome[r.chrom][r.start : r.end]
        except Exception as e:  # noqa: BLE001
            raise ValueError(f"Failed fetching {r.fasta_header} from genome FASTA {genome_fasta}") from e
        seqs.append(str(seq))
    return seqs


def load_sequences_from_fasta(sequences_fasta, records):  # type: (Union[str, Path], Sequence[EnhancerRecord]) -> List[str]
    """
    Load sequences from a FASTA whose headers match `chrom:start-end` for each record.

    This is the most portable mode: users can generate this FASTA with any tool they like
    and then let BRAIN-MAGNET handle splitting + target writing.
    """
    sequences_fasta = Path(sequences_fasta)
    needed = set([r.fasta_header for r in records])
    found = {}  # type: dict

    header = None
    seq_chunks = []  # type: List[str]
    with sequences_fasta.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # flush previous
                if header is not None:
                    if header in needed and header not in found:
                        found[header] = "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
                continue
            seq_chunks.append(line)

        # flush last
        if header is not None:
            if header in needed and header not in found:
                found[header] = "".join(seq_chunks)

    missing = [h for h in needed if h not in found]
    if missing:
        raise ValueError(
            "Missing {n} sequences in {p}. Example missing header: {ex}".format(
                n=len(missing), p=sequences_fasta, ex=missing[0]
            )
        )

    return [found[r.fasta_header] for r in records]


def split_indices(n, seed, train_frac, valid_frac):  # type: (int, int, float, float) -> Tuple[List[int], List[int], List[int]]
    if n <= 0:
        return [], [], []
    if not (0 < train_frac < 1) or not (0 < valid_frac < 1) or (train_frac + valid_frac) >= 1:
        raise ValueError("train_frac and valid_frac must be in (0,1) and sum to < 1")

    import random

    idx = list(range(n))
    rng = random.Random(seed)
    rng.shuffle(idx)

    n_train = int(n * train_frac)  # floor; matches the observed historical split sizes
    n_valid = int(n * valid_frac)  # floor
    train_idx = idx[:n_train]
    valid_idx = idx[n_train : n_train + n_valid]
    test_idx = idx[n_train + n_valid :]
    return train_idx, valid_idx, test_idx


T = TypeVar("T")


def subset_by_index(items, indices):  # type: (Sequence[T], Sequence[int]) -> List[T]
    return [items[i] for i in indices]


def prepare_enhancer_dataset(
    *,
    enhancer_activity_tsv,
    genome_fasta=None,
    sequences_fasta=None,
    out_dir,
    seed: int = 13,
    train_frac: float = 0.8,
    valid_frac: float = 0.1,
    include_reversed: bool = True,
    write_full_fasta: bool = True,
    target_columns=None,  # type: Optional[Sequence[str]]
) -> None:
    """
    Prepare FASTA + target files starting from `Enhancer_activity.txt`.

    Outputs (inside `out_dir`):
    - Enhancer.fa (optional; forward only)
    - train_set/
        - Sequences_Train.fa
        - Sequences_Valid.fa
        - Sequences_Test.fa
        - Sequences_activity_Train.txt
        - Sequences_activity_Valid.txt
        - Sequences_activity_Test.txt

    The split is reproducible via `seed`. Within each split FASTA, the historical format is:
    - forward entries first
    - then the reverse-complement entries appended, with headers like `>START-END Reversed:`
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    records, activity_names = read_enhancer_activity_tsv(enhancer_activity_tsv)
    if sequences_fasta is not None:
        try:
            # Preferred: `sequences_fasta` contains one record per interval,
            # with headers like `>chr10:89391824-89392824`.
            sequences = load_sequences_from_fasta(sequences_fasta, records)
        except ValueError:
            # Convenience: if someone accidentally passes a *genome* FASTA instead,
            # try interval extraction from the same FASTA (requires pyfaidx).
            sequences = fetch_genome_sequences(sequences_fasta, records)
    elif genome_fasta is not None:
        sequences = fetch_genome_sequences(genome_fasta, records)
    else:
        raise ValueError("Provide either `sequences_fasta` or `genome_fasta` to obtain sequences.")

    if write_full_fasta:
        write_fasta(iter_fasta(records, sequences, include_reversed=False), out_dir / "Enhancer.fa")

    train_idx, valid_idx, test_idx = split_indices(len(records), seed, train_frac, valid_frac)

    train_recs = subset_by_index(records, train_idx)
    valid_recs = subset_by_index(records, valid_idx)
    test_recs = subset_by_index(records, test_idx)

    train_seqs = subset_by_index(sequences, train_idx)
    valid_seqs = subset_by_index(sequences, valid_idx)
    test_seqs = subset_by_index(sequences, test_idx)

    train_dir = out_dir / "train_set"
    write_fasta(
        iter_fasta(train_recs, train_seqs, include_reversed=include_reversed),
        train_dir / "Sequences_Train.fa",
    )
    write_fasta(
        iter_fasta(valid_recs, valid_seqs, include_reversed=include_reversed),
        train_dir / "Sequences_Valid.fa",
    )
    write_fasta(
        iter_fasta(test_recs, test_seqs, include_reversed=include_reversed),
        train_dir / "Sequences_Test.fa",
    )

    if target_columns is None:
        activity_indices = list(range(len(activity_names)))
    else:
        name_to_idx = {n: i for i, n in enumerate(activity_names)}
        missing = [c for c in target_columns if c not in name_to_idx]
        if missing:
            raise ValueError("Requested target columns not found in file header: {m}".format(m=missing))
        activity_indices = [name_to_idx[c] for c in target_columns]

    write_activity_tsv(
        train_recs,
        train_dir / "Sequences_activity_Train.txt",
        include_reversed=include_reversed,
        activity_names=activity_names,
        activity_indices=activity_indices,
    )
    write_activity_tsv(
        valid_recs,
        train_dir / "Sequences_activity_Valid.txt",
        include_reversed=include_reversed,
        activity_names=activity_names,
        activity_indices=activity_indices,
    )
    write_activity_tsv(
        test_recs,
        train_dir / "Sequences_activity_Test.txt",
        include_reversed=include_reversed,
        activity_names=activity_names,
        activity_indices=activity_indices,
    )

