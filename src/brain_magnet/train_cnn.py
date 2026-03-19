import csv
import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


def _ensure_torch():
    _ensure_ld_library_path()
    try:
        import torch  # noqa: F401
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(
            "PyTorch is required for training.\n"
            "Install it in your environment (e.g. pixi/conda/pip) and rerun.\n"
            "Original import error: {e}".format(e=e)
        )


def _ensure_ld_library_path():
    prefix = os.environ.get("CONDA_PREFIX") or os.environ.get("VIRTUAL_ENV") or sys.prefix
    candidates = [os.path.join(prefix, "lib"), os.path.join(prefix, "lib64")]
    libs = [p for p in candidates if os.path.isdir(p)]
    if not libs:
        return

    cur = os.environ.get("LD_LIBRARY_PATH", "")
    parts = [p for p in cur.split(":") if p]
    for lib in libs:
        if lib not in parts:
            parts.insert(0, lib)
    os.environ["LD_LIBRARY_PATH"] = ":".join(parts)


def read_fasta_sequences(path: Path, uppercase: bool = True) -> Tuple[List[str], List[str]]:
    headers = []
    seqs = []
    h = None
    chunks = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    s = "".join(chunks)
                    if uppercase:
                        s = s.upper()
                    headers.append(h)
                    seqs.append(s)
                h = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if h is not None:
            s = "".join(chunks)
            if uppercase:
                s = s.upper()
            headers.append(h)
            seqs.append(s)
    return headers, seqs


def one_hot_encode(sequence: str, alphabet: str = "ACGT", neutral: str = "N") -> "np.ndarray":
    import numpy as np
    # A,C,G,T -> one-hot; N (and any other character) -> zeros row
    table = np.zeros((256, 4), dtype=np.float32)
    for i, ch in enumerate(alphabet):
        table[ord(ch)][i] = 1.0
    arr = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
    return table[arr]


def pad_to_length(x: "np.ndarray", target_len: int) -> "np.ndarray":
    import numpy as np
    # x: [L, 4]
    L = x.shape[0]
    if L == target_len:
        return x
    if L > target_len:
        start = (L - target_len) // 2
        return x[start : start + target_len]
    diff = target_len - L
    left = diff // 2
    right = diff - left
    return np.pad(x, [(left, right), (0, 0)], mode="constant")


def load_activity_column(path: Path, target_column: str) -> "np.ndarray":
    import numpy as np
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Empty header in {p}".format(p=path))
        if target_column not in reader.fieldnames:
            raise ValueError(
                "Target column '{c}' not found in {p}. Available: {cols}".format(
                    c=target_column, p=path, cols=reader.fieldnames
                )
            )
        values = []
        for row in reader:
            values.append(float(row[target_column]))
    return np.asarray(values, dtype=np.float32)


def filter_forward_only(
    headers: Sequence[str], seqs: Sequence[str], y: "np.ndarray"
) -> Tuple[List[str], List[str], "np.ndarray"]:
    import numpy as np
    keep_h = []
    keep_s = []
    keep_y = []
    for h, s, v in zip(headers, seqs, y):
        if "Reversed" in h:
            continue
        keep_h.append(h)
        keep_s.append(s)
        keep_y.append(v)
    return keep_h, keep_s, np.asarray(keep_y, dtype=np.float32)


def make_X(seqs: Sequence[str], seq_len: int) -> "np.ndarray":
    import numpy as np
    X = np.zeros((len(seqs), seq_len, 4), dtype=np.float32)
    for i, s in enumerate(seqs):
        X[i] = pad_to_length(one_hot_encode(s), seq_len)
    # replace any weird numeric issues
    return np.nan_to_num(X)


def pearsonr_np(x: "np.ndarray", y: "np.ndarray") -> float:
    import numpy as np
    x = x.astype(np.float64)
    y = y.astype(np.float64)
    x = x - x.mean()
    y = y - y.mean()
    denom = np.sqrt((x * x).sum()) * np.sqrt((y * y).sum())
    if denom == 0:
        return float("nan")
    return float((x * y).sum() / denom)


def spearmanr_np(x: "np.ndarray", y: "np.ndarray") -> float:
    import numpy as np
    # rank with average for ties
    def rankdata(a):
        temp = a.argsort()
        ranks = np.empty_like(temp, dtype=np.float64)
        ranks[temp] = np.arange(len(a), dtype=np.float64)
        # average ties
        _, inv, counts = np.unique(a, return_inverse=True, return_counts=True)
        sums = np.bincount(inv, ranks)
        avg = sums / counts
        return avg[inv]

    rx = rankdata(x)
    ry = rankdata(y)
    return pearsonr_np(rx, ry)


class TrainConfig:
    def __init__(
        self,
        train_set_dir: Path,
        target_column: str,
        output_dir: Path,
        seq_len: int = 1000,
        batch_size: int = 128,
        epochs: int = 100,
        patience: int = 20,
        learning_rate: float = 1e-2,
        seed: int = 913,
        device: str = "auto",
        drop_reverse: bool = True,
    ):
        self.train_set_dir = train_set_dir
        self.target_column = target_column
        self.output_dir = output_dir
        self.seq_len = seq_len
        self.batch_size = batch_size
        self.epochs = epochs
        self.patience = patience
        self.learning_rate = learning_rate
        self.seed = seed
        self.device = device
        self.drop_reverse = drop_reverse


def train_and_eval(cfg: TrainConfig) -> Dict[str, Dict[str, float]]:
    _ensure_ld_library_path()
    import numpy as np  # type: ignore
    _ensure_torch()
    import torch
    from torch import nn
    from torch.nn import BatchNorm1d, BatchNorm2d, Conv2d, Dropout, Flatten, Linear, MaxPool2d, ReLU, Sequential
    from torch.utils.data import DataLoader, TensorDataset

    # output layout
    model_dir = cfg.output_dir / "models"
    preds_dir = cfg.output_dir / "preds_targets" / cfg.target_column
    model_dir.mkdir(parents=True, exist_ok=True)
    preds_dir.mkdir(parents=True, exist_ok=True)

    # seed
    np.random.seed(cfg.seed)
    torch.manual_seed(cfg.seed)

    if cfg.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(cfg.device)

    # The CNN architecture pools the sequence down by 2^3 = 8 and then flattens
    # with 512 channels and 1 height dimension: flat_size = 512 * (seq_len // 8).
    # seq_len must be divisible by 8.
    if cfg.seq_len % 8 != 0:
        raise ValueError(
            "--seq-len must be divisible by 8 (got {}).".format(cfg.seq_len)
        )
    cnn_flat_size = 512 * (cfg.seq_len // 8)

    def load_split(split: str):
        fasta = cfg.train_set_dir / ("Sequences_" + split + ".fa")
        act = cfg.train_set_dir / ("Sequences_activity_" + split + ".txt")
        headers, seqs = read_fasta_sequences(fasta, uppercase=True)
        y = load_activity_column(act, cfg.target_column)
        if cfg.drop_reverse:
            headers, seqs, y = filter_forward_only(headers, seqs, y)
        X = make_X(seqs, cfg.seq_len)
        # torch: [N, 4, 1, L]
        tensor_x = torch.tensor(X).permute(0, 2, 1).unsqueeze(2)
        tensor_y = torch.tensor(y).unsqueeze(1)
        ds = TensorDataset(tensor_x, tensor_y)
        dl = DataLoader(ds, batch_size=cfg.batch_size, shuffle=(split == "Train"))
        return headers, seqs, y, dl

    _, _, _, train_dl = load_split("Train")
    _, _, _, valid_dl = load_split("Valid")
    test_headers, test_seqs, test_y, test_dl = load_split("Test")

    class CNN_STARR(nn.Module):
        def __init__(self):
            super(CNN_STARR, self).__init__()
            self.model = Sequential(
                Conv2d(4, 128, (1, 11), padding="same"),
                BatchNorm2d(128),
                ReLU(),
                MaxPool2d((1, 2), (1, 2)),
                Conv2d(128, 256, (1, 9), padding="same"),
                BatchNorm2d(256),
                ReLU(),
                MaxPool2d((1, 2), (1, 2)),
                Conv2d(256, 512, (1, 7), padding="same"),
                BatchNorm2d(512),
                ReLU(),
                MaxPool2d((1, 2), (1, 2)),
                Flatten(),
                Linear(cnn_flat_size, 1024),
                BatchNorm1d(1024),
                ReLU(),
                Dropout(0.4),
                Linear(1024, 1024),
                BatchNorm1d(1024),
                ReLU(),
                Dropout(0.4),
                Linear(1024, 1),
            )

        def forward(self, x):
            return self.model(x)

    model = CNN_STARR().to(device)
    loss_fn = torch.nn.MSELoss().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=cfg.learning_rate)

    class EarlyStopping:
        def __init__(self, patience: int, verbose: bool, task: str):
            self.patience = patience
            self.verbose = verbose
            self.task = task
            self.counter = 0
            self.best_score = None
            self.early_stop = False
            self.val_loss_min = float("inf")

        def __call__(self, val_loss: float, model):
            score = -val_loss
            if self.best_score is None:
                self.best_score = score
                self.save_checkpoint(val_loss, model)
                return
            if score < self.best_score:
                self.counter += 1
                if self.verbose:
                    print("EarlyStopping counter: {} out of {}".format(self.counter, self.patience))
                if self.counter >= self.patience:
                    self.early_stop = True
                return
            self.best_score = score
            self.save_checkpoint(val_loss, model)
            self.counter = 0

        def save_checkpoint(self, val_loss: float, model):
            if self.verbose:
                print(
                    "Validation loss decreased ({:.5f} --> {:.5f}). Saving model ...".format(
                        self.val_loss_min, val_loss
                    )
                )
            ckpt = model_dir / "checkpoint_{}.pth".format(self.task)
            torch.save(model.state_dict(), str(ckpt))
            self.val_loss_min = val_loss

    early = EarlyStopping(patience=cfg.patience, verbose=True, task=cfg.target_column)

    history = {"train_loss": [], "valid_loss": []}
    for epoch in range(cfg.epochs):
        model.train()
        train_losses = []
        for inputs, targets in train_dl:
            inputs = inputs.to(device)
            targets = targets.to(device)
            outputs = model(inputs)
            loss = loss_fn(outputs, targets)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            train_losses.append(float(loss.item()))

        model.eval()
        valid_losses = []
        with torch.no_grad():
            for inputs, targets in valid_dl:
                inputs = inputs.to(device)
                targets = targets.to(device)
                outputs = model(inputs)
                loss = loss_fn(outputs, targets)
                valid_losses.append(float(loss.item()))

        train_loss = float(np.mean(train_losses)) if train_losses else float("nan")
        valid_loss = float(np.mean(valid_losses)) if valid_losses else float("nan")
        history["train_loss"].append(train_loss)
        history["valid_loss"].append(valid_loss)

        print("[{} / {}] train_loss {:.6f} valid_loss {:.6f}".format(epoch + 1, cfg.epochs, train_loss, valid_loss))
        early(valid_loss, model)
        if early.early_stop:
            print("Early stopping")
            break

    # load best checkpoint
    ckpt = model_dir / "checkpoint_{}.pth".format(cfg.target_column)
    model.load_state_dict(torch.load(str(ckpt), map_location=device, weights_only=False))

    def predict(dl) -> Tuple[np.ndarray, np.ndarray]:
        model.eval()
        preds = []
        targets = []
        with torch.no_grad():
            for inputs, y in dl:
                inputs = inputs.to(device)
                out = model(inputs).detach().cpu().numpy().squeeze(1)
                preds.append(out)
                targets.append(y.detach().cpu().numpy().squeeze(1))
        return np.concatenate(targets), np.concatenate(preds)

    metrics = {}
    test_targets_arr = None
    test_preds_arr = None
    for split, dl in [("Train", train_dl), ("Valid", valid_dl), ("Test", test_dl)]:
        t, p = predict(dl)
        mse = float(np.mean((t - p) ** 2))
        pr = pearsonr_np(t, p)
        sr = spearmanr_np(t, p)
        metrics[split] = {"mse": mse, "pearson": pr, "spearman": sr}

        np.save(str(preds_dir / "targets_{}_{}.npy".format(cfg.target_column, split)), t)
        np.save(str(preds_dir / "preds_{}_{}.npy".format(cfg.target_column, split)), p)

        if split == "Test":
            test_targets_arr = t
            test_preds_arr = p

    # save merged test table (forward-only) for convenience
    out_csv = preds_dir / "test_{}.tsv".format(cfg.target_column)
    with out_csv.open("w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["location", "sequence", "target", "pred"])
        for h, s, tval, pval in zip(test_headers, test_seqs, test_targets_arr, test_preds_arr):
            w.writerow([h, s, float(tval), float(pval)])

    # save metrics + history
    (cfg.output_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))
    (cfg.output_dir / "history.json").write_text(json.dumps(history, indent=2))
    print("Saved:", cfg.output_dir / "metrics.json")
    print("Saved:", cfg.output_dir / "history.json")
    print("Saved:", out_csv)

    return metrics

