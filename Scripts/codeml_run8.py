import os
from gwf import Workflow
from gwf import AnonymousTarget

gwf = Workflow(defaults={'cores':2,'memory':'32gb','account':'genomics_diapause','walltime':'96:00:00'})


# -------------------------
# Prep target genes for codeml using codeml_prep7.py
# -------------------------
def codeml_prep_target(gene, gene_dir, tree, templates_dir):
    description = f"Prep target {gene} for codeml"
    #define input paths
    site_temp = os.path.join(templates_dir, "site_template.txt")  
    phylo311 = "/home/anenaur/miniforge3/envs/phylo311/bin/python"

    #define output paths
    phy = os.path.join(gene_dir, f"{gene}.phy")
    nwk = os.path.join(gene_dir, f"{gene}.nwk")
    site_ctl = os.path.join(gene_dir, f"{gene}_site.ctl")
    log = os.path.join(gene_dir, f"{gene}_prep.log")

    #define workflow
    inputs = [tree, site_temp]
    outputs = [phy, nwk, site_ctl, log]
    options = {}
    spec = f"""
    {phylo311} codeml_prep7.py \
    --gene {gene} \
    --gene_dir {gene_dir} \
    --tree {tree} \
    --site_template {site_temp}
    > {log} 2>&1
    """
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
# Summarize results using summarize6.py
# -------------------------
def summarize(gene, gene_dir, result_dir):
    description = f"Make summary of codeml output for {gene}"
    
    #define input and output paths
    site = os.path.join(gene_dir, f"{gene}_site.mlc")
    out = os.path.join(result_dir, f"{gene}_summary.txt")
    python = "/home/anenaur/miniforge3/envs/codeml/bin/python"

    #define workflow
    inputs = [site]
    outputs = [out]
    options = {}
    spec = f"""
    {python} summarize6.py \
    --gene {gene} \
    --site {site} \
    --out {out} \
    --threshold 5.99
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options={},spec=spec)


# -------------------------
# Extract beb using extract_beb2.py
# -------------------------
def extract_beb(gene, gene_dir, result_dir):
    description = f"Extract BEB from codeml output alternative model for {gene} if lnL is significant"
    
    #define input and output paths
    site = os.path.join(gene_dir, f"{gene}_site.mlc")
    out = os.path.join(result_dir, f"{gene}_beb.txt")
    python = "/home/anenaur/miniforge3/envs/codeml/bin/python"

    #define workflow
    inputs = [site]
    outputs = [out]
    options = {}
    spec = f"""
    {python} extract_beb2.py \
    --site {site} \
    --out {out} \
    --threshold 5.99 \
    --name {gene}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options={},spec=spec)


#Define paths
base = os.getcwd()
data = os.path.join(base, "8data")
result = os.path.join(base, "8results")
tree = os.path.join(base, "full_omm_tree.nwk")
gene_list = ["IGF1", "COX5A", "QRICH1", "STC1", "MYO1A", "DAP", "IGFBP5", "CAPZB", "STC2"]

#Execute workflow on test genes
for gene in gene_list:
    #define path for directory with input data
    gene_dir = os.path.join(data, gene)
    
    #execute preperation script
    gwf.target_from_template(
        f"prep_{gene}_8",
        codeml_prep_target(gene, gene_dir, tree, base)) 
    
    #run codeml on site models   
    gwf.target_from_template(
        f"codeml_site_{gene}_8",
        run_codeml(gene, gene_dir, "site"))
    
    #make summary and extract beb from codeml output
    gwf.target_from_template(
        f"summarize_{gene}_8",
        summarize(gene, gene_dir, result))
    gwf.target_from_template(
        f"exctract_beb_{gene}_8",
        extract_beb(gene, gene_dir, result)
    )
