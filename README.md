Project Overview

This project investigates the genome-wide effects of histone deacetylase inhibition (HDACi) on histone acetylation in the mouse hippocampus. Using ChIP-Seq, we mapped four histone marks:

H3K4me3
AcH3K9,14
AcH4K12
pan-AcH2B

Mice were treated with either Trichostatin A (TSA, 2.4 mg/kg) or Vehicle (DMSO/Saline), and hippocampal tissue was collected 30 minutes post-treatment.
The goal is to identify promoter regions with hyperacetylation induced by HDAC inhibition and assess the impact on gene regulation.

Data Sources:

Raw and processed ChIP-Seq data are included:
GSE43439_RAW.tar – Raw WIG files
GSE43439_gbd1_3_X_SICER-W200-G600-E10-islandfiltered.wig.normalized.wig.gz – Normalized and filtered WIG data

Mouse genome annotation:

GTF file: Mus_musculus.GRCm39.109.gtf (for promoters)

Analysis Overview:

Read WIG files for TSA and Vehicle samples for all four histone marks.
Define promoter regions based on mouse GTF annotation.
Calculate average ChIP signal per promoter region.
Compute Fold Change (TSA vs Vehicle) for each histone mark.
Visualize results using heatmaps to identify hyperacetylated and hypoacetylated genes.
