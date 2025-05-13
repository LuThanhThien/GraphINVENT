import pandas as pd

def get_sub_data(path, out, max_line=500_000):
    new_lines = []
    with open(path, "r") as f:
        print("Reading file...")
        # read the first line to get the header
        # yeild f.readline()
        header = f.readline()
        # read the rest of the file to get the number of lines
        num_line = 0
        for line in f:
            num_line += 1
            if num_line % 1000 == 0:
                print("Reading line: ", num_line)
            if num_line > max_line + 1:
                break
            new_lines.append(line)

    with open(out, "w") as f:
        print("Writing file...")
        # write the header
        f.write(header)
        # write the rest of the file
        for line in new_lines:
            f.write(line)
    print("Done!")

def tsv2csv(path, out):
    file_tsv = pd.read_csv(path, sep="\t")
    file_tsv.to_csv(out, sep=",", index=False)
    
def generate_smi_file(path, out):
    # Read in the CSV
    df = pd.read_csv(path, sep=",", dtype=str)

    # Keep only the SMILES column as a DataFrame
    df = df[["SMILES"]]

    # Reset the integer index into a column, rename it "Name"
    df = df.reset_index().rename(columns={"index": "Name"})

    # Reorder so SMILES is first and Name second (reset_index put Name first by default)
    df = df[["SMILES", "Name"]]

    # Write out space-delimited, no extra index
    df.to_csv(out, sep=" ", index=False)
    print("Done!")

# pre process data:
# 
# get_sub_data("data/Papyrus/05.5_combined_set_without_stereochemistry.tsv", 
#         "data/Papyrus/05.5_combined_set_without_stereochemistry_mini.tsv")

# tsv2csv("data/Papyrus/05.5_combined_set_without_stereochemistry_mini.tsv", 
#         "data/Papyrus/05.5_combined_set_without_stereochemistry_mini.csv")

generate_smi_file("data/Papyrus/05.5_combined_set_without_stereochemistry_mini.csv",
                  "data/Papyrus/05.5_combined_set_without_stereochemistry_mini.smi")

