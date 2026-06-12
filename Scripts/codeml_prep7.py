from Bio import AlignIO
from ete3 import Tree
from Bio import Phylo
from Bio import SeqIO
from Bio.AlignIO import PhylipIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from io import StringIO
import re
import os
import argparse

# -------------------------
# FASTA to phylip
# -------------------------
def fasta_to_phylip_codeml(gene, wrk_dir, kept_species):
    base = os.getcwd() #save current directory
    os.chdir(wrk_dir) #change working directory

    #define in- and output files
    fasta_file = f"{gene}.fasta"
    phylip_file = f"{gene}.phy"

    #read sequences
    records = list(SeqIO.parse(fasta_file, "fasta")) #get species name - making sure name isnt truncated as standard phylip only allows 10 character names
    original_count = len(records) #count number of species before filtering

    #remove species not in subset
    records = [r for r in records if r.id in kept_species]

    #remove sequences with >75% gaps - otherwise codeml crashes
    records = [r for r in records if (str(r.seq).count("-") / len(r.seq)) <= 0.75]
    print(f"Removed {original_count - len(records)} species from {gene} due to >75% gaps (or not in selected species)")#print how many were removed

    num_seqs = len(records)
    seq_len = len(records[0].seq)

    #check sequences same length, otherwise codeml throws an error
    if not all(len(r.seq) == seq_len for r in records):
        raise ValueError(f"Sequences in {fasta_file} must be same length")

    #ID width: longest ID + 2 spaces
    max_id_len = max(len(r.id) for r in records)
    id_width = max_id_len + 2 #ensuring all aligned sequences start at the exact same place, while ensuring at least two spaces between species name and sequence, as codeml otherwise dont get that the name is longer than 10 characters

    #write phylip file
    with open(phylip_file, "w") as out:
        #phylip header
        out.write(f"{num_seqs} {seq_len}\n") #number of species + alignment length

        #write sequences
        for r in records:
            #pad species names to same length, add sequence immediately after
            padded_id = r.id.ljust(id_width)
            out.write(f"{padded_id}{r.seq}\n")

    os.chdir(base) #return to original directory
    print(f"PHYLIP file for {gene} written to {phylip_file}")

# -------------------------
# Unroot tree, prune to match alignment
# -------------------------
def newick_unroot_with_foreground(newick_tree, gene, data_dir):
    #define in- and output paths
    tree = Phylo.read(newick_tree, "newick")
    output_tree = f"{data_dir}/{gene}/{gene}.nwk"
    phylip_file = f"{data_dir}/{gene}/{gene}.phy"  # filtered alignment

    #get species from filtered alignment
    alignment = AlignIO.read(phylip_file, "phylip-relaxed")
    kept_species = {record.id for record in alignment}

    #create file like object in memory without saving to disk for easy processing
    handle = StringIO() 
    Phylo.write(tree, handle, "newick") 

    #load/create tree object with ete3
    t = Tree(handle.getvalue())

    #find all species in tree
    tree_species = {leaf.name for leaf in t.iter_leaves()}
    
    #find which species are not in the alignment
    to_remove = tree_species - kept_species

    if to_remove:
        print(f"{gene}: removing {len(to_remove)} species from tree")

    #prune tree (keeps only species in alignment)
    t.prune(kept_species, preserve_branch_length=False)

    #unroot tree
    t.unroot()

    #write tree to newick file
    t.write(format=9, outfile=output_tree) #format tells what whould be included in the tree, 9 means only leaf names and no branch lengths

    print(f"{gene}: tree pruned to {len(t)} species and written to {output_tree}")


# -------------------------
# Control files
# -------------------------
def make_codeml_ctl(template_ctl, new_ctl, seqfile, outfile, treefile):
    #define paths to be added to the control file
    replacements = {
        "seqfile": seqfile,
        "outfile": outfile,
        "treefile": treefile
    }

    #create copy of template file and add defined paths
    with open(template_ctl) as infile, open(new_ctl, "w") as out:
        for line in infile:
            for key, value in replacements.items():
                if value and re.match(rf"\s*{key}\s*=", line):
                    line = f"{key} = {value}\n"
            out.write(line)


# -------------------------
# Prepare site model control file
# -------------------------
def make_site_ctl(gene, template, data_dir):
    #path for sequence file
    seqfile = os.path.join(data_dir, f"{gene}.phy")

    #paths for new control files
    ctl1 = os.path.join(data_dir, f"{gene}_site.ctl")

    #output files for codeml
    out1 = os.path.join(data_dir, f"{gene}_site.mlc") 
    
    #paths for tree file
    tree = os.path.join(data_dir, f"{gene}.nwk")

    #make control file for site model
    make_codeml_ctl(
        template,
        ctl1,
        seqfile,
        out1,
        tree
    )


#Make script executable in gwf
if __name__ == "__main__":
    import argparse

    #define required arguments for execution
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", required=True)
    parser.add_argument("--gene_dir", required=True)
    parser.add_argument("--tree", required=True)
    parser.add_argument("--site_template", required=True)
   

    #read the arguments provided
    args = parser.parse_args()
    gene = args.gene
    gene_dir = args.gene_dir
    tree = args.tree
    site_temp = args.site_template

    #list of kept species
    species = [ "Choloepus_didactylus", "Dasypus_novemcinctus", "Trichechus_manatus_latirostris",
               "Elephas_maximus_indicus", "Elephantulus_edwardii", "Echinops_telfairi", "Meriones_unguiculatus", 
               "Rattus_norvegicus", "Rattus_rattus", "Grammomys_surdaster", "Mastomys_coucha", "Mus_pahari",
               "Mus_musculus", "Mus_caroli", "Onychomys_torridus", "Peromyscus_californicus_insignis",
               "Peromyscus_maniculatus_bairdii", "Peromyscus_leucopus", "Mesocricetus_auratus", "Arvicola_amphibius",
               "Microtus_fortis", "Talpa_occidentalis", "Condylura_cristata", "Erinaceus_europaeus", "Suncus_etruscus",
               "Sorex_araneus", "Canis_lupus_dingo", "Vulpes_vulpes", "Meles_meles", "Lontra_canadensis",
               "Enhydra_lutris_kenyoni", "Neogale_vison", "Mustela_erminea", "Ailuropoda_melanoleuca", "Ursus_americanus",
               "Ursus_maritimus", "Ursus_arctos", "Odobenus_rosmarus_divergens", "Callorhinus_ursinus", 
               "Zalophus_californianus", "Eumetopias_jubatus", "Halichoerus_grypus", "Neomonachus_schauinslandi", 
               "Leptonychotes_weddellii", "Mirounga_leonina", "Mirounga_angustirostris"]

    #execute function converting alignment fasta to phylip
    fasta_to_phylip_codeml(
        gene, 
        gene_dir,
        species
    )

    #execute function filtering and unrooting tree
    newick_unroot_with_foreground(
        tree,
        gene,
        os.path.dirname(gene_dir)
    )

    #execute function making codeml control file for site model
    make_site_ctl(
        gene,
        site_temp,
        gene_dir
    )

