import re
import argparse

# -------------------------
# Read file
# -------------------------
def read_file(file):
    with open(file) as f:
        return f.read()


# -------------------------
# Extract general info
# -------------------------
def parse_general_info(file):
    #read input file
    text = read_file(file)

    #find number of species (ns), length of alignment (ls), and codon model in input file
    ns_match = re.search(r"ns\s*=\s*(\d+)", text)
    ls_match = re.search(r"ls\s*=\s*(\d+)", text)
    codon_match = re.search(r"Codon frequency model:\s*(.+)", text)

    n_species = ns_match.group(1) if ns_match else "NA"
    aln_length = ls_match.group(1) if ls_match else "NA"
    codon_freq = codon_match.group(1).strip() if codon_match else "NA"

    return n_species, aln_length, codon_freq


# -------------------------
# Extract lnL
# -------------------------
def get_lnL(file):
    #read input file
    text = read_file(file)

    #find lnL in input file
    match = re.search(r"lnL\(.*?\):\s*([-\d\.]+)", text)
    if match:
        return float(match.group(1))

    return float("nan")  # prevents crashing


# -------------------------
# Extract kappa
# -------------------------
def get_kappa(file):
    #read input file
    text = read_file(file)

    #find kappa in input file
    match = re.search(r"kappa\s*\(ts/tv\)\s*=\s*([-\d\.]+)", text)
    return match.group(1) if match else "NA"


# -------------------------
# Extract omega
# -------------------------
def get_omega_table(file):
    #read input file
    text = read_file(file)

    # ---- Case 1: Branch model ----
    #find omega in input file - values following "w (dN/dS) for branches"
    branch_match = re.search(
        r"w\s*\(dN/dS\)\s*for branches:\s*([\d\.\s]+)",
        text
    )

    if branch_match:
        values = branch_match.group(1).strip().split()
        if len(values) >= 3:
            return [
                f"background_w\t{values[0]}",
                f"foreground1_w\t{values[1]}",
                f"foreground2_w\t{values[2]}"
            ]
    
    # ---- Case 2: single omega (M0 / null model) ----
    #find omega value in input file - following "omega (dN/dS)"
    single_match = re.search(
        r"omega\s*\(dN/dS\)\s*=\s*([\d\.]+)",
        text
    )

    if single_match:
        return [f"omega\t{single_match.group(1)}"]

    # ---- Case 3: fallback ----
    else:
        return ["no omega found"]


# -------------------------
# Write model info
# -------------------------
def write_model_info(out, file, model_name):
    out.write(f"---{model_name}---\n")

    kappa = get_kappa(file)
    out.write(f"kappa\t{kappa}\n")

    omega_table = get_omega_table(file)
    for line in omega_table:
        out.write(line + "\n")

    lnL = get_lnL(file)
    out.write(f"lnL\t{lnL:.5f}\n\n")

    return lnL


# -------------------------
# Main summary function
# -------------------------
def summarize(gene, M0_file, branch_file, out_file, threshold=5.99):
    #default threshold of 5.99 as 2Deltal should be compared to chi squared with 2 degrees of freedom
    n_species, aln_length, codon_freq = parse_general_info(branch_file)

    with open(out_file, "w") as out:

        #general info
        out.write(f"Gene\t{gene}\n")
        out.write(f"Number of Species\t{n_species}\n")
        out.write(f"Alignment length (Number of codons)\t{aln_length}\n")
        out.write(f"Codon model\t{codon_freq}\n\n")

        #branch test
        lnL_branch = write_model_info(out, branch_file, "Branch alternative")
        lnL_M0 = write_model_info(out, M0_file, "Branch null")

        out.write("*branch model test*\n")
        delta2 = 2 * (lnL_branch - lnL_M0)
        out.write(f"2delta_lnL\t{delta2:.5f}\n")

        if delta2 >= threshold:
            out.write("SIGNIFICANT\n")
        else:
            out.write("NON-SIGNIFICANT\n")


# -------------------------
# CLI
# -------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    #define required arguments for execution
    parser.add_argument("--gene", required=True)
    parser.add_argument("--M0", required=True)
    parser.add_argument("--branch", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--threshold", type=float, default=5.99)

    #read the arguments provided
    args = parser.parse_args()
    gene = args.gene
    M0 = args.M0
    branch = args.branch
    out = args.out
    threshold = args.threshold

    #execute summary function with provided arguments
    summarize(
        gene,
        M0,
        branch,
        out,
        threshold
    )