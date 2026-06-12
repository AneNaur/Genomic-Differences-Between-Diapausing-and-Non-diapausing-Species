import os
from gwf import Workflow
from gwf import AnonymousTarget

gwf = Workflow(defaults={'cores':2,'memory':'16gb','account':'genomics_diapause','walltime':'96:00:00'})


# -------------------------
# Prep target genes for codeml using codeml_prep5.py
# -------------------------
def codeml_prep_target(gene, gene_dir, tree, templates_dir):
    description = f"Prep target {gene} for codeml"
    #define input paths
    branch_temp = os.path.join(templates_dir, "branch_template.txt")
    M0_temp = os.path.join(templates_dir, "M0_template.txt")    
    phylo311 = "/home/anenaur/miniforge3/envs/phylo311/bin/python"

    #define output paths
    phy = os.path.join(gene_dir, f"{gene}.phy")
    nwk = os.path.join(gene_dir, f"{gene}.nwk")
    branch_ctl = os.path.join(gene_dir, f"{gene}_branch.ctl")
    M0_ctl = os.path.join(gene_dir, f"{gene}_M0.ctl")
    log = os.path.join(gene_dir, f"{gene}_prep.log")

    #define workflow
    inputs = [tree, branch_temp, M0_temp]
    outputs = [phy, nwk, branch_ctl, M0_ctl, log]
    options = {}
    spec = f"""
    {phylo311} codeml_prep5.py \
    --gene {gene} \
    --gene_dir {gene_dir} \
    --tree {tree} \
    --branch_template {branch_temp} \
    --M0_template {M0_temp} \
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
    outputs = [out, log]  # Include log as output
    options = {}
    spec = f"""
    cd {gene_dir}
    codeml {gene}_{model}.ctl | tee {gene}_{model}.log
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# -------------------------
# Summarize results using summarize4.py
# -------------------------
def summarize(gene, gene_dir, result_dir):
    description = f"Make summary of codeml output for {gene}"
    #define input and output paths
    M0 = os.path.join(gene_dir, f"{gene}_M0.mlc")
    branch = os.path.join(gene_dir, f"{gene}_branch.mlc") 
    out = os.path.join(result_dir, f"{gene}_summary.txt")
    python = "/home/anenaur/miniforge3/envs/codeml/bin/python"

    #define workflow
    inputs = [M0, branch]
    outputs = [out]
    options = {}
    spec = f"""
    {python} summarize4.py \
    --gene {gene} \
    --M0 {M0} \
    --branch {branch} \
    --out {out} \
    --threshold 3.84
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options={},spec=spec)


#Define paths
base = os.getcwd()
data = os.path.join(base, "test_data2") #"6data")
result = os.path.join(base, "test_results2") #"6results")
tree = os.path.join(base, "full_omm_tree.nwk")
gene_list1 = ["IGF1", "IGF2", "IGF1R", "IGF2R", "IGFBP2", "IGFBP4", "IGFBP5", "PAPPA", "PAPPA2", "STC1", "STC2", "NANOG", "POU5F1", "SOX2"]
gene_list2 = ["HDHD2", "HADH", "HK2", "PSMA7", "COX17", "COX5A", "ABAT", "DAP", "CAPZB", "PMPCA", "ACTR1A", "DDIT4", "DDOST", "EHF", "EVPL", "GAS7", "MAP1B", "MRPL27", "MYO1A", "NPAS4", "PAFAH1B3", "PAQR3", "PHF24", "POLR1C", "QRICH1", "TGFB3", "TGM1", "USH1C", "ZC3HC1", "ZCCHC17"]

#Execute workflow
for gene in gene_list1:
    #define path for directory with input data
    gene_dir = os.path.join(data, gene)
    
    #execute preperation script
    gwf.target_from_template(
        f"prep_{gene}_6",
        codeml_prep_target(gene, gene_dir, tree, base)) 
    
    #run codeml on alternative and null model for branch model test
    gwf.target_from_template(
        f"codeml_branch_{gene}_6",
        run_codeml(gene, gene_dir, "branch"))    
    gwf.target_from_template(
        f"codeml_M0_{gene}_6",
        run_codeml(gene, gene_dir, "M0"))
     
    #summarize results
    gwf.target_from_template(
        f"summarize_{gene}_6",
        summarize(gene, gene_dir, result))


#Execute workflow
for gene in gene_list2:
    #define path for directory with input data
    gene_dir = os.path.join(data, gene)
    
    #execute preperation script
    gwf.target_from_template(
        f"prep_{gene}_6",
        codeml_prep_target(gene, gene_dir, tree, base)) 
    
    #run codeml on alternative and null model for branch model test
    gwf.target_from_template(
        f"codeml_branch_{gene}_6",
        run_codeml(gene, gene_dir, "branch"))    
    gwf.target_from_template(
        f"codeml_M0_{gene}_6",
        run_codeml(gene, gene_dir, "M0"))
     
    #summarize results
    gwf.target_from_template(
        f"summarize_{gene}_6",
        summarize(gene, gene_dir, result))

