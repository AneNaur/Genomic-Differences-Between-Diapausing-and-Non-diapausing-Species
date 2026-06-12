import re
import argparse
import math
import os

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
    return match.group(1) if match else "kappa\tNA"


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
        if len(values) >= 2:
            return [
                f"background_w\t{values[0]}",
                f"foreground_w\t{values[1]}"
            ]

    # ---- Case 2: Branch-site models ----
    lines = text.split("\n")
    omega_lines = []
    capture = False

    for line in lines:
        #find start of omega table in input file - following "MLEs of dN/dS (w) for site classes (K=4)"
        if re.search(r"MLEs of dN/dS \(w\)", line):
            omega_lines.append(line.strip())
            capture = True
            continue

        if capture:
            stripped = line.strip()

            #find relevant lines: "site class", "proportion", "background w", and "foreground w"
            if re.match(r"(site class|proportion|background\s+w|foreground\s+w)", stripped):
                omega_lines.append(stripped)

            #stop after foreground row
            if stripped.startswith("foreground"):
                break

    if omega_lines:
        cleaned_lines = []

        for line in omega_lines:
            line = line.strip()

            #fix labels BEFORE tabbing
            line = re.sub(r"background\s+w", "background_w", line)
            line = re.sub(r"foreground\s+w", "foreground_w", line)
            line = re.sub(r"site\s+class", "site_class", line)

            #insert tabs between each table entry
            line = re.sub(r"\s+", "\t", line)

            cleaned_lines.append(line)

        return cleaned_lines

    # ---- Case 3: single omega (M0 / null model) ----
    #find omega value in input file - following "omega (dN/dS)"
    single_match = re.search(
        r"omega\s*\(dN/dS\)\s*=\s*([\d\.]+)",
        text
    )

    if single_match:
        return [f"omega\t{single_match.group(1)}"]

    # ---- Case 4: fallback ----
    else:
        return ["no omega found"]


# -------------------------
# Write model info
# -------------------------
def write_model_info(out, file, model_name):
    out.write(f"---{model_name}---\n")

    if not os.path.exists(file):
        out.write("FILE NOT FOUND\n\n")
        return float("nan")

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
def summarize(gene, null_file, alt_file, M0_file, branch_file, out_file, threshold=3.84):

    n_species, aln_length, codon_freq = parse_general_info(alt_file)

    with open(out_file, "w") as out:

        #general info
        out.write(f"Gene\t{gene}\n")
        out.write(f"Number of Species\t{n_species}\n")
        out.write(f"Alignment length (Number of codons)\t{aln_length}\n")
        out.write(f"Codon model\t{codon_freq}\n\n")

        #branch-site test
        lnL_alt = write_model_info(out, alt_file, "Branch-site alternative")
        lnL_null = write_model_info(out, null_file, "Branch-site null")

        out.write("*branch-site model test*\n")
        if not (math.isnan(lnL_alt) or math.isnan(lnL_null)):
            delta1 = 2 * (lnL_alt - lnL_null)
            out.write(f"2delta_lnL\t{delta1:.5f}\n")

            if delta1 >= threshold:
                out.write("SIGNIFICANT\n\n")
            else:
                out.write("NON-SIGNIFICANT\n\n")
        else:
            out.write("missing lnL\n\n")

        #branch test
        lnL_branch = write_model_info(out, branch_file, "Branch alternative")
        lnL_M0 = write_model_info(out, M0_file, "Branch null")

        out.write("*branch model test*\n")
        if not (math.isnan(lnL_branch) or math.isnan(lnL_M0)):
            delta2 = 2 * (lnL_branch - lnL_M0)
            out.write(f"2delta_lnL\t{delta2:.5f}\n")
            if delta2 >= threshold:
                out.write("SIGNIFICANT\n")
            else:
                out.write("NON-SIGNIFICANT\n")            
        else: 
            out.write("missing lnL\n")



# -------------------------
# CLI
# -------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    #define required arguments for execution
    parser.add_argument("--gene", required=True)
    parser.add_argument("--null", required=True)
    parser.add_argument("--alt", required=True)
    parser.add_argument("--M0", required=True)
    parser.add_argument("--branch", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--threshold", type=float, default=3.84)

    #read the arguments provided
    args = parser.parse_args()
    gene = args.gene
    null = args.null
    alt = args.alt
    M0 = args.M0
    branch = args.branch
    out = args.out
    threshold = args.threshold

    #execute summary function with provided arguments
    summarize(
        gene,
        null,
        alt,
        M0,
        branch,
        out,
        threshold
    )