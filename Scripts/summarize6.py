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
# Split input file into model chunks
# -------------------------
def split_model_sections(file):

    # read entire file
    text = read_file(file)

    # split before each NSsites model header
    sections = re.split(
        r'(?=^\s*NSsites Model\s+\d+:)',
        text,
        flags=re.MULTILINE
    )

    # remove empty chunks
    sections = [s.strip() for s in sections if s.strip()]

    return sections


# -------------------------
# Extract model
# -------------------------
def get_models(text):
    match = re.search(r'NSsites\s+(Model\s+\d+:\s*[^(\n]+)', text)
    if match:
        return match.group(1).strip()
    
    return "NA"


# -------------------------
# Extract kappa
# -------------------------
def get_kappa(text):
    #find kappa in input chunck
    match = re.search(r"kappa\s*\(ts/tv\)\s*=\s*([-\d\.]+)", text)
    if match:
        return float(match.group(1))
    
    return float("nan")


# -------------------------
# Extract number of classes K
# -------------------------
def get_K(line):
    match = re.search(r"K\s*=\s*(\d+)", line)
    return match.group(1) if match else "?"


# -------------------------
# Extract omega
# -------------------------
def get_omega(text):
    lines = text.split("\n")
    omega_lines = []
    capture = False
    k = None

    for line in lines:
        if re.search(r"MLEs of dN/dS \(w\)", line):
            k = get_K(line)
            omega_lines.append(f"Omega for site classes (K={k})")
            capture = True
            continue

        if capture:
            stripped = line.strip()

            #find relevant lines: "p:" and "w:"
            if re.match(r"^(p|w)\s*:", stripped):
                omega_lines.append(stripped)

            #stop after foreground row
            if stripped.startswith("w:"):
                break
    
    return omega_lines


# -------------------------
# Extract lnL
# -------------------------
def get_lnL(text):
    #find lnL in input file
    match = re.search(r"lnL\(.*?\):\s*([-\d\.]+)", text)
    if match:
        return float(match.group(1))

    return float("nan")  # prevents crashing


# -------------------------
# Main summary function
# -------------------------
def summarize(gene, in_file, out_file, threshold=5.99):
    #default threshold of 5.99 as 2Deltal should be compared to chi squared with 2 degrees of freedom
    n_species, aln_length, codon_freq = parse_general_info(in_file)

    with open(out_file, "w") as out:
        #general info
        out.write(f"Gene\t{gene}\n")
        out.write(f"Number of Species\t{n_species}\n")
        out.write(f"Alignment length (Number of codons)\t{aln_length}\n")
        out.write(f"Codon model\t{codon_freq}\n\n")

        chunks = split_model_sections(in_file)
        m = len(chunks)
        lnLs = []

        for i in range(1,m):
            #model name
            chunk = chunks[i]
            model_name = get_models(chunk)
            out.write(f"---{model_name}---\n")
            #kappa
            kappa = get_kappa(chunk)
            out.write(f"kappa\t{kappa}\n")
            #omega
            omega = get_omega(chunk)
            for line in omega:
                out.write(line + "\n")
            #lnL
            lnL = get_lnL(chunk)
            out.write(f"lnL\t{lnL}\n\n")
            lnLs.append(lnL)
        
        #2deltalnL
        out.write("*site model test: M1a vs. M2a*\n")
        delta1 = 2 * (lnLs[1] - lnLs[0])
        out.write(f"2delta_lnL\t{delta1:.5f}\n")
        if delta1 >= threshold:
                out.write("SIGNIFICANT\n\n")
        else:
            out.write("NON-SIGNIFICANT\n\n")
        out.write("*site model test: M7 vs. M8*\n")
        delta2 = 2 * (lnLs[3] - lnLs[2])
        out.write(f"2delta_lnL\t{delta2:.5f}\n")
        if delta2 >= threshold:
                out.write("SIGNIFICANT")
        else:
            out.write("NON-SIGNIFICANT")        


# -------------------------
# CLI
# -------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    #define required arguments for execution
    parser.add_argument("--gene", required=True)
    parser.add_argument("--site", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--threshold", type=float, default=5.99)

    #read the arguments provided
    args = parser.parse_args()
    gene = args.gene
    site = args.site
    out = args.out
    threshold = args.threshold

    #execute summary function with provided arguments
    summarize(
        gene,
        site,
        out,
        threshold
    )