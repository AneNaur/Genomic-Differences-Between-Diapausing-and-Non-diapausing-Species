import os
from gwf import Workflow
from gwf import AnonymousTarget
import json

gwf = Workflow(defaults={'cores':2,'memory':'32gb','account':'genomics_diapause','walltime':'72:00:00'})


# -------------------------
# Prep target genes for codeml using codeml_prep4.py
# -------------------------
def codeml_prep_target(gene, gene_dir, tree, templates_dir, order):
    description = f"Prep target {gene} for codeml"
    #define input paths
    alt_temp = os.path.join(templates_dir, "alt_template.txt")
    null_temp = os.path.join(templates_dir, "null_template.txt")
    branch_temp = os.path.join(templates_dir, "branch_template.txt")
    M0_temp = os.path.join(templates_dir, "M0_template.txt")    
    phylo311 = "/home/anenaur/miniforge3/envs/phylo311/bin/python"

    #define output paths
    phy = os.path.join(gene_dir, f"{gene}_{order}.phy")
    nwk = os.path.join(gene_dir, f"{gene}_{order}.nwk")
    alt_ctl = os.path.join(gene_dir, f"{gene}_{order}_alt.ctl")
    null_ctl = os.path.join(gene_dir, f"{gene}_{order}_null.ctl")
    branch_ctl = os.path.join(gene_dir, f"{gene}_{order}_branch.ctl")
    M0_ctl = os.path.join(gene_dir, f"{gene}_{order}_M0.ctl")
    log = os.path.join(gene_dir, f"{gene}_{order}_prep.log")

    #define workflow
    inputs = [tree, alt_temp, null_temp, branch_temp, M0_temp]
    outputs = [phy, nwk, alt_ctl, null_ctl, branch_ctl, M0_ctl, log]
    options = {}
    spec = f"""
    {phylo311} codeml_prep4.py \
    --gene {gene} \
    --gene_dir {gene_dir} \
    --tree {tree} \
    --alt_template {alt_temp} \
    --null_template {null_temp} \
    --branch_template {branch_temp} \
    --M0_template {M0_temp} \
    --order {order}
    > {log} 2>&1
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# -------------------------
# Run codeml + make log file
# -------------------------
def run_codeml(gene, gene_dir, model, order):
    description = f"Run codeml ({model}) for {gene} in {order}"
    
    #define input and output paths
    ctl = os.path.join(gene_dir, f"{gene}_{order}_{model}.ctl")
    out = ctl.replace(".ctl", ".mlc")
    log = os.path.join(gene_dir, f"{gene}_{order}_{model}.log")

    #define workflow
    inputs = [ctl]
    outputs = [out, log]
    options = {}
    spec = f"""
    cd {gene_dir}
    codeml {gene}_{order}_{model}.ctl | tee {gene}_{order}_{model}.log
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# -------------------------
# Summarize results using summarize3.py
# -------------------------
def summarize(gene, order_list, gene_dir, result_dir):
    description = f"Make summary of codeml output for {gene}"
    #define input and output paths
    null_list = {}
    alt_list = {}
    M0_list = {}
    branch_list = {}
    out = os.path.join(result_dir, f"{gene}_summary.txt")
    python = "/home/anenaur/miniforge3/envs/codeml/bin/python"    
    
    for order in order_list:
        #define input and output paths
        null_list[order] = os.path.join(gene_dir, f"{gene}_{order}_null.mlc")
        alt_list[order] = os.path.join(gene_dir, f"{gene}_{order}_alt.mlc")
        M0_list[order] = os.path.join(gene_dir, f"{gene}_{order}_M0.mlc")
        branch_list[order] = os.path.join(gene_dir, f"{gene}_{order}_branch.mlc")

    #define workflow
    inputs = [null_list, alt_list, M0_list, branch_list]
    outputs = [out]
    options = {}
    spec = f"""
    {python} summarize3.py \
    --gene {gene} \
    --order_list '{json.dumps(order_list)}' \
    --null_list '{json.dumps(null_list)}' \
    --alt_list '{json.dumps(alt_list)}' \
    --M0_list '{json.dumps(M0_list)}' \
    --branch_list '{json.dumps(branch_list)}' \
    --out {out} \
    --threshold 3.84
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options={},spec=spec)


# -------------------------
# Extract beb using extract_beb1.py
# -------------------------
def extract_beb(gene, gene_dir, order, result_dir):
    description = f"Extract BEB from codeml output alternative model for {gene} if lnL is significant"
    
    #define input and output paths
    null = os.path.join(gene_dir, f"{gene}_{order}_null.mlc")
    alt = os.path.join(gene_dir, f"{gene}_{order}_alt.mlc")
    out = os.path.join(result_dir, f"{gene}_{order}_beb.txt")
    python = "/home/anenaur/miniforge3/envs/codeml/bin/python"

    #define workflow
    inputs = [null, alt]
    outputs = [out]
    options = {}
    spec = f"""
    {python} extract_beb1.py \
    --gene {gene}_{order} \
    --null {null} \
    --alt {alt} \
    --out {out} \
    --threshold 3.84 \
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options={},spec=spec)


#Define paths
base = os.getcwd()
data = os.path.join(base, "4data")
result = os.path.join(base, "4results")
tree = os.path.join(base, "full_omm_tree.nwk")
gene_list = ["IGF1", "IGF2", "IGF1R", "IGF2R", "IGFBP2", "IGFBP4", "IGFBP5", "PAPPA", "PAPPA2", "STC1", "STC2"]


#Execute workflow
for gene in gene_list2:
    #define path for directory with input data
    gene_dir = os.path.join(data, gene)
    for order in order_list:
        #execute preperation script
        gwf.target_from_template(
            f"prep_{gene}_{order}_4",
            codeml_prep_target(gene, gene_dir, tree, base, order)) 
    
        #run codeml on all four models
        gwf.target_from_template(
            f"codeml_alt_{gene}_{order}_4",
            run_codeml(gene, gene_dir, "alt", order))
        gwf.target_from_template(
            f"codeml_null_{gene}_{order}_4",
            run_codeml(gene, gene_dir, "null", order))
        gwf.target_from_template(
            f"codeml_branch_{gene}_{order}_4",
            run_codeml(gene, gene_dir, "branch", order))    
        gwf.target_from_template(
            f"codeml_M0_{gene}_{order}_4",
            run_codeml(gene, gene_dir, "M0", order))
        
        #extract beb from codeml output
        gwf.target_from_template(
            f"exctract_beb_{gene}_{order}_4",
            extract_beb(gene, gene_dir, order, result))        

    #make summary of results
    gwf.target_from_template(
        f"summarize_{gene}_4",
        summarize(gene, order_list, gene_dir, result))
    
