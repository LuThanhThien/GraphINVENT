from pathlib import Path
import os
import argparse


def main():
    p = argparse.ArgumentParser(
        description="Split a .smi file into train, val, and test sets."
    )
    p.add_argument("input",  help=".smi files")
    p.add_argument(
        "--train",
        type=float,
        default=0.8,
        help="fraction of data to use for training (default: 0.8)",
    )
    p.add_argument(
        "--test",
        type=float,
        default=0.1,
        help="fraction of data to use for testing (default: 0.1)",
    )
    p.add_argument(
        "--output",
        "-o",
        type=str,
        default="data/",
        help="output directory (default: current directory)",
    )
    args = p.parse_args()
    input_file = Path(args.input)
    output_dir = Path(args.output)
    
    if not output_dir.exists():
        os.makedirs(output_dir)
    if not input_file.exists():
        print(f"File {input_file} does not exist.")
        return
    if not input_file.is_file():
        print(f"{input_file} is not a file.")
        return
    
    with open(input_file) as fin:
        lines = fin.readlines()
    
    new_lines = set()
    for i, line in enumerate(lines):
        if line.startswith("SMILES"):
            continue
        new_lines.add(line.strip().split()[0] + " 1\n")
    
    new_lines = list(new_lines)
    num_lines = len(new_lines)
    train_size = int(num_lines * args.train)
    test_size = int(num_lines * args.test)
    val_size = num_lines - train_size - test_size
    train_lines = new_lines[:train_size]
    test_lines = new_lines[train_size:train_size + test_size]
    val_lines = new_lines[train_size + test_size:]
    
    with open(output_dir / "train.smi", "w") as fout:
        fout.writelines(train_lines)
    with open(output_dir / "test.smi", "w") as fout:
        fout.writelines(test_lines)
    with open(output_dir / "val.smi", "w") as fout:
        fout.writelines(val_lines)
        
    print(f"Split {num_lines} lines into {len(train_lines)} train, {len(test_lines)} test, and {len(val_lines)} val lines.")

if __name__ == "__main__":
    main()