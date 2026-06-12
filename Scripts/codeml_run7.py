import os
from gwf import Workflow
from gwf import AnonymousTarget

gwf = Workflow(defaults={'cores':2,'memory':'32gb','account':'genomics_diapause','walltime':'96:00:00'})


# -------------------------
# Prep target genes for codeml using codeml_prep6.py
# -------------------------
def codeml_prep_target(gene, gene_dir, tree, templates_dir, foreground):
    description = f"Prep target {gene} for codeml"
    #define input paths
    alt_template = os.path.join(templates_dir, "alt_template.txt")
    null_template = os.path.join(templates_dir, "null_template.txt")
    phylo311 = "/home/anenaur/miniforge3/envs/phylo311/bin/python"

    #define output paths
    phy = os.path.join(gene_dir, f"{gene}.phy")
    nwk = os.path.join(gene_dir, f"{gene}.nwk")
    alt_ctl = os.path.join(gene_dir, f"{gene}_alt.ctl")
    null_ctl = os.path.join(gene_dir, f"{gene}_null.ctl")
    log = os.path.join(gene_dir, f"{gene}_prep.log")

    #define workflow
    inputs = [tree, alt_template, null_template]
    outputs = [phy, nwk, alt_ctl, null_ctl, log]
    options = {}
    spec = f"""
    {phylo311} codeml_prep6.py \
    --gene {gene} \
    --gene_dir {gene_dir} \
    --tree {tree} \
    --alt_template {alt_template} \
    --null_template {null_template} \
    --foreground {foreground} \
    > {log} 2>&1
    #"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# -------------------------
# Run codeml + make log file
# -------------------------
def run_codeml(gene, gene_dir, model):
    description = f"Run codeml ({model}) for {gene}"
    
    #define input and output paths
    ctl = os.path.join(gene_dir, f"{gene}_{model}.ctl")
    out = ctl.replace(".ctl", ".mlc")
    log = os.path.join(gene_dir, f"{gene}_{model}.log")

    #define workflow
    inputs = [ctl]
    outputs = [out, log]
    options = {}
    spec = f"""
    cd {gene_dir}
    codeml {gene}_{model}.ctl | tee {gene}_{model}.log
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# -------------------------
# Summarize results using summarize5.py
# -------------------------
def summarize(gene, gene_dir, result_dir, foreground):
    description = f"Make summary of codeml output for {gene}"
    
    #define input and output paths
    null = os.path.join(gene_dir, f"{gene}_null.mlc")
    alt = os.path.join(gene_dir, f"{gene}_alt.mlc")
    out = os.path.join(result_dir, f"{gene}_summary.txt")
    python = "/home/anenaur/miniforge3/envs/codeml/bin/python"

    #define workflow
    inputs = [null, alt]
    outputs = [out]
    options = {}
    spec = f"""
    {python} summarize5.py \
    --gene {gene} \
    --null {null} \
    --alt {alt} \
    --out {out} \
    --foreground {foreground}
    --threshold 3.84
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options={},spec=spec)

# -------------------------
# Extract beb using extract_beb1.py
# -------------------------
def extract_beb(gene, gene_dir, result_dir):
    description = f"Extract BEB from codeml output alternative model for {gene} if lnL is significant"
    
    #define input and output paths
    null = os.path.join(gene_dir, f"{gene}_null.mlc")
    alt = os.path.join(gene_dir, f"{gene}_alt.mlc")
    out = os.path.join(result_dir, f"{gene}_beb.txt")
    python = "/home/anenaur/miniforge3/envs/codeml/bin/python"

    #define workflow
    inputs = [null, alt]
    outputs = [out]
    options = {}
    spec = f"""
    {python} extract_beb1.py \
    --null {null} \
    --alt {alt} \
    --out {out} \
    --threshold 3.84
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options={},spec=spec)


#Define paths
base = os.getcwd()
data = os.path.join(base, "7data")
result = os.path.join(base, "7results")
tree = os.path.join(base, "full_omm_tree.nwk")
alt_temp = os.path.join(base, "alt_template.txt")
null_temp = os.path.join(base, "null_template.txt")
gene_list1 = ["SOX2_1", "NANOG_1", "POU5F1_1", "IGF1", "IGFBP5", "STC2"] #Obligate
gene_list2 = ["SOX2_2", "NANOG_2", "POU5F1_2", "IGF1R", "IGF2R", "IGFBP4", "PAPPA2", "STC1"] #Facultative


#Execute workflow
for gene in gene_list1:
    #define path for directory with input data
    gene_dir = os.path.join(data, gene)
    
    #execute preperation script
    gwf.target_from_template(
        f"prep_{gene}_7",
        codeml_prep_target(gene, gene_dir, tree, base, "Obligate")) 
    
    #run codeml on alternative and null model
    gwf.target_from_template(
        f"codeml_alt_{gene}_7",
        run_codeml(gene, gene_dir, "alt"))
    gwf.target_from_template(
        f"codeml_null_{gene}_7",
        run_codeml(gene, gene_dir, "null"))

    #make summary and extract beb from codeml output
    gwf.target_from_template(
        f"summarize_{gene}_7",
        summarize(gene, gene_dir, result, "Obligate"))
    gwf.target_from_template(
        f"exctract_beb_{gene}_1",
        extract_beb(gene, gene_dir, result)
    )

for gene in gene_list2:
#define path for directory with input data
    gene_dir = os.path.join(data, gene)
    
    #execute preperation script
    gwf.target_from_template(
        f"prep_{gene}_7",
        codeml_prep_target(gene, gene_dir, tree, base, "Facultative")) 
    
    #run codeml on alternative and null model
    gwf.target_from_template(
        f"codeml_alt_{gene}_7",
        run_codeml(gene, gene_dir, "alt"))
    gwf.target_from_template(
        f"codeml_null_{gene}_7",
        run_codeml(gene, gene_dir, "null"))

    #make summary and extract beb from codeml output
    gwf.target_from_template(
        f"summarize_{gene}_7",
        summarize(gene, gene_dir, result, "Facultative"))
    gwf.target_from_template(
        f"exctract_beb_{gene}_1",
        extract_beb(gene, gene_dir, result))    

