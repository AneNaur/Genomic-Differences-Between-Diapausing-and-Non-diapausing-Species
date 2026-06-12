import re
import argparse
import json


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
    text = read_file(file)

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
    text = read_file(file)
    match = re.search(r"lnL\(.*?\):\s*([-\d\.]+)", text)

    return float(match.group(1)) if match else float("nan")


# -------------------------
# Extract kappa
# -------------------------
def get_kappa(file):
    text = read_file(file)
    match = re.search(r"kappa\s*\(ts/tv\)\s*=\s*([-\d\.]+)", text)

    return match.group(1) if match else "NA"


# -------------------------
# Extract omega
# -------------------------
def get_omega_table(file):
    text = read_file(file)

    # Case 1: branch model
    branch_match = re.search(
        r"w\s*\(dN/dS\)\s*for branches:\s*([\d\.\s]+)",
        text
    )

    if branch_match:
        values = branch_match.group(1).split()
        if len(values) >= 2:
            return [
                f"background_w\t{values[0]}",
                f"foreground_w\t{values[1]}"
            ]

    # Case 2: branch-site model
    lines = text.split("\n")
    omega_lines = []
    capture = False

    for line in lines:
        if re.search(r"MLEs of dN/dS \(w\)", line):
            capture = True
            continue

        if capture:
            stripped = line.strip()

            if re.match(r"(site class|proportion|background\s+w|foreground\s+w)", stripped):
                omega_lines.append(stripped)

            if stripped.startswith("foreground"):
                break

    if omega_lines:
        cleaned = []
        for line in omega_lines:
            line = re.sub(r"background\s+w", "background_w", line)
            line = re.sub(r"foreground\s+w", "foreground_w", line)
            line = re.sub(r"site\s+class", "site_class", line)
            line = re.sub(r"\s+", "\t", line)
            cleaned.append(line)
        return cleaned

    # Case 3: single omega
    single = re.search(r"omega\s*\(dN/dS\)\s*=\s*([\d\.]+)", text)
    if single:
        return [f"omega\t{single.group(1)}"]

    # Case 4: fallback
    else:
        return ["no omega found"]


# -------------------------
# Write model block
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
# One order summary
# -------------------------
def summarize_order(out, order, null_file, alt_file, M0_file, branch_file, threshold):
    out.write(f"\n========== {order} ==========\n")

    # branch-site test
    lnL_alt = write_model_info(out, alt_file, "Branch-site alternative")
    lnL_null = write_model_info(out, null_file, "Branch-site null")

    delta1 = 2 * (lnL_alt - lnL_null)
    out.write("*branch-site model test*\n")
    out.write(f"2delta_lnL\t{delta1:.5f}\n")
    out.write("SIGNIFICANT\n\n" if delta1 >= threshold else "NON-SIGNIFICANT\n\n")

    # branch test
    lnL_branch = write_model_info(out, branch_file, "Branch alternative")
    lnL_M0 = write_model_info(out, M0_file, "Branch null")

    delta2 = 2 * (lnL_branch - lnL_M0)
    out.write("*branch model test*\n")
    out.write(f"2delta_lnL\t{delta2:.5f}\n")
    out.write("SIGNIFICANT\n\n" if delta2 >= threshold else "NON-SIGNIFICANT\n\n")


# -------------------------
# Main summary
# -------------------------
def summarize2(gene, order_list, null_list, alt_list, M0_list, branch_list, out_file, threshold):

    main_file = alt_list[order_list[0]]
    n_species, aln_length, codon_freq = parse_general_info(main_file)

    with open(out_file, "w") as out:
        out.write(f"Gene\t{gene}\n")
        out.write(f"Alignment length\t{aln_length}\n")
        out.write(f"Codon model\t{codon_freq}\n")

        for order in order_list:
            summarize_order(
                out,
                order,
                null_list[order],
                alt_list[order],
                M0_list[order],
                branch_list[order],
                threshold
            )


# -------------------------
# CLI
# -------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--gene", required=True)
    parser.add_argument("--order_list", required=True)

    parser.add_argument("--null_list", required=True)
    parser.add_argument("--alt_list", required=True)
    parser.add_argument("--M0_list", required=True)
    parser.add_argument("--branch_list", required=True)

    parser.add_argument("--out", required=True)
    parser.add_argument("--threshold", type=float, default=3.84)

    args = parser.parse_args()

    # JSON parsing (IMPORTANT FIX)
    order_list = json.loads(args.order_list)
    null_list = json.loads(args.null_list)
    alt_list = json.loads(args.alt_list)
    M0_list = json.loads(args.M0_list)
    branch_list = json.loads(args.branch_list)

    summarize2(
        args.gene,
        order_list,
        null_list,
        alt_list,
        M0_list,
        branch_list,
        args.out,
        args.threshold
    )