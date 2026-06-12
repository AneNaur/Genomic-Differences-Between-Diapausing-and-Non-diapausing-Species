import re
import argparse


# -------------------------
# Read file
# -------------------------
def read_file(file):
    with open(file) as f:
        return f.read()


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
# Extract BEB
# -------------------------

def get_beb(text):
    lines = text.splitlines()

    beb_start = None
    beb = []

    for i, line in enumerate(lines):
        if "Bayes Empirical Bayes" in line:
            beb_start = i
            break

    if beb_start is None:
        return []

    for line in lines[beb_start:]:
        if "The grid" in line:
            break

        if re.match(r"\s*\d+\s+[A-Z-]\s+\d+\.\d+", line):
            beb.append(line.strip())

    return beb

#def get_beb(text):
    lines = text.splitlines()

    beb_start = None
    beb = []

    for i, line in enumerate(lines):
        if "Bayes Empirical Bayes" in line:
            beb_start = i
            break

    if beb_start is None:
        return []

    for line in lines[beb_start:]:
        if "The grid" in line:
            break

        if re.match(r"\s*\d+\s+[A-Z-]\s+\d+\.\d+", line):
            beb.append(line.strip())

    return beb

#def get_beb(text):
    lines = text.splitlines()

    beb_start = None
    beb = []

    for i, line in enumerate(lines):
        if "Bayes Empirical Bayes" in line:
            beb_start = i
            break

    if beb_start is None:
        return []

    for line in lines[beb_start:]:
        if "The grid" in line:
            break

        if re.match(r"\s*\d+\s+[A-Z]\s+", line):
            beb.append(line.strip())

    return beb


def summarize_beb(gene, null_file, alt_file, out_file, threshold):
    null_text = read_file(null_file)
    alt_text = read_file(alt_file)

    lnL_null = get_lnL(null_text)
    lnL_alt = get_lnL(alt_text)
    delta = 2 * (lnL_alt - lnL_null) #calculate LRT statistic

    #make output file
    with open(out_file, "w") as out:
        out.write(f"---{gene}---\n")
        out.write(f"Null lnL\t{lnL_null}\n")
        out.write(f"Alternative lnL\t{lnL_alt}\n")
        out.write(f"2delta_lnL\t{delta:.5f}\n") #write LRT statistic

        #check if significant
        if delta >= threshold:
            out.write("SIGNIFICANT\n")
            out.write("site\tAA\tprob\n")

            for site in get_beb(alt_text):
                fields = site.split()
                out.write("\t".join(fields) + "\n")
        
        else:
            out.write("NON-SIGNIFICANT\n")
            out.write("No beb to report")


        
        

        #read codeml output file for alternative model
        
        #with open(alt_file) as f:
        #    lines = f.readlines()

        #beb_start = None
        #for i, line in enumerate(lines):
        #    if "Bayes Empirical Bayes" in line:
        #        beb_start = i
        #        break

        #if beb_start is None:
        #    out.write("NO_BEB_FOUND\n")
        #    return

        #for line in lines[beb_start:]:
        #    if "The grid" in line:
        #        break
        #    if re.match(r"\s*\d+\s+[A-Z]\s+", line):
        #        out.write(line.strip() + "\n")


#Make script executable in gwf
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    #define required arguments for execution
    parser.add_argument("--gene", required=True)
    parser.add_argument("--null", required=True)
    parser.add_argument("--alt", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--threshold", type=float, default=3.84) #if not provided default to significance level 5%

    #read the arguments provided
    args = parser.parse_args()
    gene = args.gene
    null = args.null
    alt = args.alt
    out = args.out
    threshold = args.threshold

    #execute function extracting beb values
    summarize_beb(gene, null, alt, out, threshold)