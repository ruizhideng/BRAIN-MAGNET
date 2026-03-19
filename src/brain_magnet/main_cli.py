import argparse


def main(argv=None):
    p = argparse.ArgumentParser(prog="brain-magnet", description="BRAIN-MAGNET command line tools.")
    sub = p.add_subparsers(dest="command")

    # -------------------------
    # prepare_data subcommand
    # -------------------------
    p_prep = sub.add_parser(
        "prepare_data",
        help="Prepare Enhancer.fa + train/valid/test splits from an activity table.",
    )
    p_prep.add_argument("--enhancer-activity", required=True, help="Path to activity TSV (BED + activity columns).")
    p_prep.add_argument(
        "--sequences-fasta",
        default=None,
        help="FASTA with headers `chrom:start-end` matching the activity table rows (recommended).",
    )
    p_prep.add_argument(
        "--genome-fasta",
        default=None,
        help="Reference genome FASTA (optional; requires pyfaidx in your env).",
    )
    p_prep.add_argument("--out-dir", required=True, help="Output directory.")
    p_prep.add_argument("--seed", type=int, default=13)
    p_prep.add_argument("--train-frac", type=float, default=0.8)
    p_prep.add_argument("--valid-frac", type=float, default=0.1)
    p_prep.add_argument("--target-cols", default=None, help="Comma-separated activity columns to write (default: all).")
    p_prep.add_argument("--no-reversed", action="store_true")
    p_prep.add_argument("--no-full-fasta", action="store_true")

    # -------------
    # train command
    # -------------
    p_train = sub.add_parser("train", help="Train the CNN model from a prepared train_set/ folder.")
    p_train.add_argument("--train-set-dir", required=True, help="Path to train_set/ directory.")
    p_train.add_argument("--target-column", required=True, help="Column name in Sequences_activity_*.txt to predict.")
    p_train.add_argument("--output-dir", required=True, help="Output directory for checkpoints/metrics/preds.")
    p_train.add_argument("--seq-len", type=int, default=1000)
    p_train.add_argument("--batch-size", type=int, default=128)
    p_train.add_argument("--epochs", type=int, default=100)
    p_train.add_argument("--patience", type=int, default=20)
    p_train.add_argument("--learning-rate", type=float, default=1e-2)
    p_train.add_argument("--seed", type=int, default=913)
    p_train.add_argument("--device", default="auto", help="auto|cpu|cuda")
    p_train.add_argument("--keep-reverse", action="store_true")

    args = p.parse_args(argv)

    if args.command is None:
        p.print_help()
        return 2

    if args.command == "prepare_data":
        # lazy import
        from pathlib import Path
        from brain_magnet.data_prep import prepare_enhancer_dataset

        def _parse_target_cols(s):
            if s is None:
                return None
            s = s.strip()
            if not s:
                return None
            return [x.strip() for x in s.split(",") if x.strip()]

        if args.sequences_fasta is None and args.genome_fasta is None:
            p_prep.error("Provide either --sequences-fasta or --genome-fasta")

        prepare_enhancer_dataset(
            enhancer_activity_tsv=Path(args.enhancer_activity),
            genome_fasta=Path(args.genome_fasta) if args.genome_fasta else None,
            sequences_fasta=Path(args.sequences_fasta) if args.sequences_fasta else None,
            out_dir=Path(args.out_dir),
            seed=args.seed,
            train_frac=args.train_frac,
            valid_frac=args.valid_frac,
            include_reversed=not args.no_reversed,
            write_full_fasta=not args.no_full_fasta,
            target_columns=_parse_target_cols(args.target_cols),
        )
        return 0

    if args.command == "train":
        from pathlib import Path
        from brain_magnet.train_cnn import TrainConfig, train_and_eval

        cfg = TrainConfig(
            train_set_dir=Path(args.train_set_dir),
            target_column=args.target_column,
            output_dir=Path(args.output_dir),
            seq_len=args.seq_len,
            batch_size=args.batch_size,
            epochs=args.epochs,
            patience=args.patience,
            learning_rate=args.learning_rate,
            seed=args.seed,
            device=args.device,
            drop_reverse=not args.keep_reverse,
        )
        train_and_eval(cfg)
        return 0

    p.error("Unknown command: {c}".format(c=args.command))

