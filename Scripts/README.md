This directory contains all scripts used for the analyses, and below is details of when each script is used. Figures were made with the rmd file "figures-diapause-thesis" using the following three inputs: omm_v12_190tax_CDS_tree.nexus (stored in directory Data), time.txt (stored in directory Results), and data_run6.txt (stored in directory Results), 

| Run | Scripts executed | Codon models | Frequency model | Species | Genes | Control template | Outputs |
|---|---|---|---|---|---|---|---|
| run1 | prep1, summarize1, beb1 | Branch-site null, Branch-site alternative | FMutSel | All | Pluripotency, control0 | alt_template, null_template | run1_summary, run1_beb |
| run2 | prep2, summarize2, beb1 | One-ratio, Fixed-ratios (2), Branch-site null, Branch-site alternative | FMutSel | Subset | Pluripotency, IGF | branch_template, M0_template, alt_template, null_template | run2_summary, run2_beb |
| run3 | prep3, summarize2, beb1 | One-ratio, Fixed-ratios (2), Branch-site null, Branch-site alternative | FMutSel | All | Pluripotency, Significant IGF run2 | branch_template, M0_template, alt_template, null_template | run3_summary, run3_beb |
| run4 | prep4, summarize3, beb1 | One-ratio, Fexed-ratios (2), Branch-site null, Branch-site alternative | FMutSel | All by order | IGF | branch_template, M0_template, alt_template, null_template | run4_summary, run4_beb |
| run5 | prep3, summarize2, beb1 | One-ratio, Fixed-ratios (2), Branch-site null, Branch-site alternative | F3x4 | Subset | IGF | branch_template2, M0_template2, alt_template2, null_template2 | run5_summary, run5_beb |
| run6 | prep5, summarize4 | One-ratio, Fixed-ratios (3) | FMutSel | All | Pluripotency, IGF, control2, Fu et al. | branch_template, M0_template | run6_summary |
| run7 | prep6, summarize5, beb1 | Branch-site alternative | FMutSel | All | Pluripotency, Significant IGF run6 | alt_template, null_template | run7_summary, run7_beb |
| run8 | prep7, summarize6, beb2 | Site models: M1a, M2a, M7, M8 | FMutSel | All | Furthest from w0=w1 line run6 | site_template | run8_summary, run8_beb |


