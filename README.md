# The transcriptome and methylome of the term human placenta in Gestational Diabetes Mellitus

### By: Marta Ibañez Lligoña

### Supervised by: Dr Amanda Sferruzzi-Perri and Dr Xiaohui Zhao.

Universitat de Barcelona (UB) and University of Cambridge (Department of Physiology, Development and Neuroscience). Sferruzzi-Perri lab.

### More information: missing information and results in this repository will be kept private until the paper is published. 


## Repository content

### Readme: 
Information about repository and folders content.

## Scripts ##
| Script file | Description | Additional information |
| ----------------------------- | ----------------------- | --------- |
| Differential_gene_expression_in-house_dataset.R  | Differential expression analysis of in-house dataset, and analysis of DEGs through enrichment analysis, and overlap with published gathered data from transcriptome/methylome-wide papers. | Needs confidential data to be run. |
| Meta_Analysis_Overlay.R | Overlay of genes found in the methylome/transcriptome-wide selected studies performed in Gestational Diabetes Mellitus Placenta and analysis of participant chatacteristics. | - |
| Omics_data_integration.R | Integration of transcriptome and methylome studies, correlation between DNA methylation and gene expression results, integration with the placental single-cell transcriptome and generation of graphs. | Needs raw data downloaded from GEO-NCBI to work. | 

## Data ##
| File | Description |
| ----------------------------- | ----------------------- |
| BMI_full_studies_Methylation.xlsx | BMI information from participants in DNA methylation selected studies. |
| BMI_full_studies_Transcriptome.xlsx | BMI information from participants in transcriptome-wide selected studies. |
| DEGs_Studies_full_list.xlsx | List of DEGs from the corresponding Gestational Diabetes Mellitus (GDM) transcriptome-wide studies. |
| DMGs_Studies_full_list.xlsx | List of DMGs corresponding to the probes identifies as differentially methylated from the corresponding GDM methylome-wide studies. |
| Transcriptome_Participant_characteristics.xlsx | Participant characteristics from transcriptome-wide selected studies. |
| Methylome_Participant_characteristics.xlsx | Participant characteristics from methylome-wide selected studies. |
| scSeq-Results.xlsx | Genes found in each cell type in the single-cell placental transcriptome. |
| Combined_Data1_Data2_matrix_transcriptome.csv | Combined normalized values from the integrated transcriptome datasets. |
| Final_Samples_Data1&2_transcriptome.csv | Table of final selected samples and their corresponding information to the integrated transcriptome-wide datasets. |

## Results ##
| File | Description |
| ----------------------------- | ----------------------- |
| Common_DMGs.xlsx | 44 DMGs found to be present in all the three methylome-wide selected studies. |
| Common_genes_in-house&published_data.xlsx | Genes found to be differentially expressed and methylated in the in-house dataset, selected transcriptome/methylome-wide studies. |
| Reactome_results_DMGs.csv | Reactome results from analysis of DMGs found to be present in all selected methylome-wide studies. |
| DEGs_data_integration.csv | DEGs from differential expression analysis of integrated transcriptome-wide datasets. |
| GO_Biological_Process_2021_table_TRANSCRIPTOME.txt | GO enrichment analysis of DEGs from the integrated selected transcriptome-wide datasets for biological process category. |
| GO_Cellular_Component_2021_table_TRANSCRIPTOME.txt | GO enrichment analysis of DEGs from the integrated selected transcriptome-wide datasets for cellular component category. |
| GO_Molecular_Function_2021_table_TRANSCRIPTOME.txt | GO enrichment analysis of DEGs from the integrated selected transcriptome-wide datasets for molecular function category. |
| KEGG_2021_Human_table_TRANSCRIPTOME.txt | KEGG pathway analysis of DEGs from the integrated selected transcriptome-wide datasets. |
| DMPs_data_integration.csv | DMPs from differential methylation analysis of integrated transcriptome-wide datasets. |
| GO_Biological_Process_2021_table_METHYLATION.txt | GO enrichment analysis of DMGs from the integrated selected methylome-wide datasets for biological process category. |
| GO_Cellular_Component_2021_table_METHYLATION.txt | GO enrichment analysis of DMGs from the integrated selected methylome-wide datasets for cellular component category. |
| GO_Molecular_Function_2021_table_METHYLATION.txt | GO enrichment analysis of DMGs from the integrated selected methylome-wide datasets for molecular function category. |
| KEGG_2021_Human_table_METHYLATION.txt | KEGG pathway analysis of DMGs from the integrated selected methylome-wide datasets. |