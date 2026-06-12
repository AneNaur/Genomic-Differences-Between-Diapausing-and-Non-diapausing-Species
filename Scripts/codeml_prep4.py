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
# FASTA to phylip with filtering based on taxonomical order
# -------------------------
def fasta_to_phylip_codeml(gene, wrk_dir, kept_species, order):
    base = os.getcwd() #save current directory
    os.chdir(wrk_dir) #change working directory

    #define in- and output files
    fasta_file = f"{gene}.fasta"
    phylip_file = f"{gene}_{order}.phy"

    #read sequences
    records = list(SeqIO.parse(fasta_file, "fasta")) #get species name - making sure name isnt truncated as standard phylip only allows 10 character names
    original_count = len(records) #count number of species before filtering

    #remove species not in the specified order
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
    print(f"PHYLIP file for {gene} and order {order} written to {phylip_file}")


# -------------------------
# Unroot tree, prune to match alignment, and tag foreground branches
# -------------------------
def newick_unroot_with_foreground(newick_tree, gene, foreground_species, data_dir, order):
    #define in- and output paths
    tree = Phylo.read(newick_tree, "newick")
    output_tree = f"{data_dir}/{gene}/{gene}_{order}.nwk"
    phylip_file = f"{data_dir}/{gene}/{gene}_{order}.phy"  # filtered alignment

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
    
    #find which species to remove
    to_remove = tree_species - kept_species

    if to_remove:
        print(f"{gene}: removing {len(to_remove)} species from tree")

    #prune tree (keeps only species in alignment)
    t.prune(kept_species, preserve_branch_length=False)

    #unroot tree
    t.unroot()

    #tag foreground species
    for leaf in t.iter_leaves():
        if leaf.name in foreground_species:
            leaf.name += " #1"

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
# Prepare branch-site model control files
# -------------------------
def make_branchsite_ctls(gene, alt_template, null_template, data_dir, order):
    #path for sequence file
    seqfile = os.path.join(data_dir, f"{gene}_{order}.phy")

    #paths for new control files
    alt_ctl = os.path.join(data_dir, f"{gene}_{order}_alt.ctl")
    null_ctl = os.path.join(data_dir, f"{gene}_{order}_null.ctl")

    #output files for codeml
    alt_out = os.path.join(data_dir, f"{gene}_{order}_alt.mlc")
    null_out = os.path.join(data_dir, f"{gene}_{order}_null.mlc")   
    
    #paths for tree file
    tree = os.path.join(data_dir, f"{gene}_{order}.nwk")

    #make control file for alternative model
    make_codeml_ctl(
        alt_template,
        alt_ctl,
        seqfile,
        alt_out,
        tree
    )

    #make control file for null model
    make_codeml_ctl(
        null_template,
        null_ctl,
        seqfile,
        null_out,
        tree
    )


# -------------------------
# Prepare branch model control files
# -------------------------
def make_branch_ctls(gene, alt_template, null_template, data_dir, order):
    #path for sequence file
    seqfile = os.path.join(data_dir, f"{gene}_{order}.phy")

    #paths for new control files
    alt_ctl = os.path.join(data_dir, f"{gene}_{order}_branch.ctl")
    null_ctl = os.path.join(data_dir, f"{gene}_{order}_M0.ctl")

    #output files for codeml
    alt_out = os.path.join(data_dir, f"{gene}_{order}_branch.mlc")
    null_out = os.path.join(data_dir, f"{gene}_{order}_M0.mlc")   
    
    #paths for tree file
    tree = os.path.join(data_dir, f"{gene}_{order}.nwk")

    #make control file for alternative model
    make_codeml_ctl(
        alt_template,
        alt_ctl,
        seqfile,
        alt_out,
        tree
    )

    #make control file for null model
    make_codeml_ctl(
        null_template,
        null_ctl,
        seqfile,
        null_out,
        tree
    )


#Make script executable in gwf
if __name__ == "__main__":
    #define required arguments for execution
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", required=True)
    parser.add_argument("--gene_dir", required=True)
    parser.add_argument("--tree", required=True)
    parser.add_argument("--alt_template", required=True)
    parser.add_argument("--null_template", required=True)
    parser.add_argument("--branch_template", required=True)
    parser.add_argument("--M0_template", required=True)
    parser.add_argument("--order", required=True)    

    #read the arguments provided
    args = parser.parse_args()
    gene = args.gene
    gene_dir = args.gene_dir
    tree = args.tree
    alt_temp = args.alt_template
    null_temp = args.null_template
    branch_temp = args.branch_template
    M0_temp = args.M0_template 
    order = args.order   

    #list of kept species
    orders = {
        "LAURASIATHERIA": ["Acinonyx_jubatus", "Ailuropoda_melanoleuca", "Artibeus_jamaicensis", "Balaenoptera_acutorostrata_scammoni",
                        "Balaenoptera_musculus", "Bison_bison_bison", "Bos_indicus", "Bos_mutus", "Bos_taurus", "Bubalus_bubalis",
                        "Budorcas_taxicolor", "Callorhinus_ursinus", "Camelus_bactrianus", "Camelus_dromedarius", "Camelus_ferus",
                        "Canis_lupus_dingo", "Capra_hircus", "Ceratotherium_simum_simum", "Cervus_canadensis", "Cervus_elaphus",
                        "Condylura_cristata", "Delphinapterus_leucas", "Desmodus_rotundus", "Enhydra_lutris_kenyoni", "Eptesicus_fuscus",
                        "Equus_asinus", "Equus_caballus", "Equus_przewalskii", "Equus_quagga", "Erinaceus_europaeus", "Eumetopias_jubatus",
                        "Felis_catus", "Globicephala_melas", "Halichoerus_grypus", "Hipposideros_armiger", "Hyaena_hyaena", "Lagenorhynchus_obliquidens",
                        "Leopardus_geoffroyi", "Leptonychotes_weddellii", "Lipotes_vexillifer", "Lontra_canadensis", "Lynx_canadensis",
                        "Manis_javanica", "Manis_pentadactyla", "Meles_meles", "Miniopterus_natalensis", "Mirounga_angustirostris",
                        "Mirounga_leonina", "Molossus_molossus", "Monodon_monoceros", "Mustela_erminea", "Mustela_putorius_furo",
                        "Myotis_brandtii", "Myotis_davidii", "Myotis_lucifugus", "Myotis_myotis", "Neogale_vison", "Neomonachus_schauinslandi",
                        "Neophocaena_asiaeorientalis_asiaeorientalis", "Odobenus_rosmarus_divergens", "Odocoileus_virginianus_texanus",
                        "Orcinus_orca", "Oryx_dammah", "Ovis_aries", "Panthera_leo", "Panthera_pardus", "Panthera_tigris", "Panthera_uncia",
                        "Phoca_vitulina", "Phocoena_sinus", "Phyllostomus_discolor", "Phyllostomus_hastatus", "Physeter_catodon",
                        "Pipistrellus_kuhlii", "Prionailurus_bengalensis", "Pteropus_alecto", "Pteropus_giganteus", "Pteropus_vampyrus",
                        "Puma_concolor", "Puma_yagouaroundi", "Rhinolophus_ferrumequinum", "Rousettus_aegyptiacus", "Sorex_araneus",
                        "Sturnira_hondurensis", "Suncus_etruscus", "Suricata_suricatta", "Sus_scrofa", "Talpa_occidentalis",
                        "Tursiops_truncatus", "Ursus_americanus", "Ursus_arctos", "Ursus_maritimus", "Vicugna_pacos",
                        "Vulpes_lagopus", "Vulpes_vulpes", "Zalophus_californianus"], 
        "GLIRES": ["Acomys_russatus", "Apodemus_sylvaticus", "Arvicanthis_niloticus", "Arvicola_amphibius",
                   "Castor_canadensis", "Cavia_porcellus", "Chinchilla_lanigera", "Cricetulus_griseus", "Dipodomys_ordii",
                   "Dipodomys_spectabilis", "Fukomys_damarensis", "Grammomys_surdaster", "Heterocephalus_glaber", 
                   "Ictidomys_tridecemlineatus", "Jaculus_jaculus", "Marmota_flaviventris", "Marmota_marmota_marmota",
                   "Marmota_monax", "Mastomys_coucha", "Meriones_unguiculatus", "Mesocricetus_auratus", "Microtus_fortis",
                   "Microtus_ochrogaster", "Microtus_oregoni", "Mus_caroli", "Mus_musculus", "Mus_pahari", "Myodes_glareolus",
                   "Nannospalax_galili", "Ochotona_curzoniae", "Ochotona_princeps", "Octodon_degus", "Onychomys_torridus",
                   "Oryctolagus_cuniculus", "Perognathus_longimembris_pacificus", "Peromyscus_californicus_insignis",
                   "Peromyscus_leucopus", "Peromyscus_maniculatus_bairdii", "Phodopus_roborovskii", "Rattus_norvegicus",
                   "Rattus_rattus", "Urocitellus_parryii"],
        "other": ["Tachyglossus_aculeatus", "Ornithorhynchus_anatinus", "Monodelphis_domestica", "Gracilinanus_agilis",
                  "Dromiciops_gliroides", "Antechinus_flavipes", "Sarcophilus_harrisii", "Trichosurus_vulpecula", "Phascolarctos_cinereus",
                  "Vombatus_ursinus", "Choloepus_didactylus", "Dasypus_novemcinctus", "Trichechus_manatus_latirostris",
                  "Elephas_maximus_indicus", "Loxodonta_africana", "Orycteropus_afer_afer", "Elephantulus_edwardii",
                  "Echinops_telfairi", "Chrysochloris_asiatica"]
                  }

    #list of foreground species
    foreground = {"Ailuropoda_melanoleuca", "Dasypus_novemcinctus",
                  "Enhydra_lutris_kenyoni", "Leptonychotes_weddellii",
                  "Meriones_unguiculatus", "Mirounga_angustirostris",
                  "Mus_musculus", "Odobenus_rosmarus_divergens",
                  "Peromyscus_maniculatus_bairdii", "Rattus_norvegicus",
                  "Sorex_araneus", "Ursus_maritimus",
                  "Zalophus_californianus", "Callorhinus_ursinus",
                  "Eumetopias_jubatus", "Halichoerus_grypus",
                  "Lontra_canadensis", "Meles_meles",
                  "Mirounga_leonina", "Mustela_erminea",
                  "Peromyscus_leucopus", "Phoca_vitulina",
                  "Ursus_americanus", "Ursus_arctos"}


    #execute function converting alignment fasta to phylip
    fasta_to_phylip_codeml(
        gene, 
        gene_dir,
        orders[order],
        order
        )

    #execute function filtering, unrooting, and tagging newick tree
    newick_unroot_with_foreground(
        tree,
        gene,
        foreground,
        os.path.dirname(gene_dir),
        order
        )

    #execute function making codeml control files for branch-site model
    make_branchsite_ctls(
        gene,
        alt_temp,
        null_temp,
        gene_dir,
        order
        )

    #execute function making codeml control files for branch model
    make_branch_ctls(
        gene,
        branch_temp,
        M0_temp,
        gene_dir,
        order
        )