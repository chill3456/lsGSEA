<p align="center">
  <a>
    <img width="40%" src="https://github.com/chill3456/lsGSEA/blob/master/assets/graphic.png" alt="graphic describing lsGSEA.R png">
  </a>
  <h1 align="center">lsGSEA</h1>
</p>




# Description 

Large scale Gene Set Enrichment Analysis (lsGSEA) uses lsGSEA.R to do Gene Set Enrichment Analysis using a variety of methods on large amounts of samples. It takes Transcripts Per Kilobase Million (TPM) gene expression counts from publically available data sets like the cBioPortal or your own TPM counts from RNA-seq or Quant-seq. It then takes the TPM counts and allows you to do ssGSEA, GSVA, zscore, or plage using gene sets from The Molecular Signatures Database (MSigDB) or your own gene sets.

# Required applications 

`install.packages(tidyverse)`
`install.packages(GSVA)`
`install.packages(dplyr)`
`install.packages(gplots)`
`install.packages(reshape2)`
`install.packages(msigdbr)`

# Input gene expression count files 

TPM gene expression counts are input of lsGSEa and can be downloaded from publically available repositories or generated from RNA-seq or Quant-seq results. 

If downloading data from the cBioPortal download files you should select "Gene Expression Quantification" "RNA-seq" and "TSV" and download the manifest. To download all samples from The Cancer Genome Atlas (TCGA) you can use this link https://portal.gdc.cancer.gov/repository?files_offset=22700&files_size=100&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22tsv%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%5B%22Gene%20Expression%20Quantification%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%5D%7D. This will give you a manifest that can be downloaded using the gdc tool.

`gdc-client download -m /Pathto/Folder/gdc_manifest.date.txt`

Each sample is output into an individual folder containing a .tsv file with gene expression counts. To combine every tsv file into a matrix that can used in lsGSEA.R you can run scripts/expression_matrix_creation_from_cbioportal.py all that is needed is to give it the path to the folder containing the cBioPortal Files. With large amounts of files this can be computationally intensive and take a while to run. 

`/scripts/expression_matrix_creation_from_cbioportal.py /example_files/cBioPortalFiles/`

It outputs a file with combined TPM counts from each downloaded sample /example_files/cBioPortalFiles/combined_gene_counts.tsv. It looks like this:

| gene_id            | gene_name | gene_type	    | 9bc8ff9e-2f36-5094-ac12-b1bf81b53c15 | 396932ab-c4b3-4526-9201-9222f5f062a5 | 53891f08-a221-4ca8-a502-cf990bdb6020 |
| ------------------ | --------- | -------------- | ------------------------------------ | ------------------------------------ | ------------------------------------ |
| ENSG00000000003.15 | TSPAN6    | protein_coding | 13.5438	                             | 40.4162                              |	0                                    |
| ENSG00000000005.6  | TNMD      | protein_coding | 0	                                   | 0.0317                               | 0                                    |
| ENSG00000000419.13 | DPM1      | protein_coding | 25.7028	                             | 105.4987	                            | 80.2825                              |
| ENSG00000000457.14 | SCYL3     | protein_coding | 2.425                                | 16.4802                              |	4.162                                |

TPMs can also be calculated from your own RNA-seq or Quant-seq. We calculate gene expression counts are caculated using featureCounts into files ending in _gene such as /example_files/making_quant_seq_TPM_files/Quant-seq-sample1_gene. Quant-seq files were generated using required_files/Homo_sapiens.GRCh37.87.conversion.txt as annotation file and can be combined into a matrix that can used in lsGSEA.R you can run scripts/getting_TPM_counts_from_Quant_seq.R all that is needed is to give it the path to the _gene Files to be used and a folder to output the results into.

`Rscript /scripts/getting_TPM_counts_from_Quant_seq.R --QuantseqGeneFiles "/example_files/quant_seq_TPM_files/Quant-seq-sample1_gene" "/example_files/quant_seq_TPM_files/Quant-seq-sample2_gene"  --outputFolder "/example_files/quant_seq_TPM_files/"`

This can also be done with RNA-seq files which were generated with required_files/gencodeV24lift37_basicbed.txt as an annotation file. It can be combined into a matrix that can used in lsGSEA.R you can run scripts/getting_TPM_counts_from_RNA_seq.R all that is needed is to give it the path to the _gene Files to be used and a folder to output the results into. 

`Rscript /scripts/getting_TPM_counts_from_RNA_seq.R --RNAseqGeneFiles "/example_files/RNA_seq_TPM_files/RNA-seq-sample1_gene" "/example_files/RNA_seq_TPM_files/RNA-seq-sample2_gene"  --outputFolder "/example_files/RNA_seq_TPM_files/"`

Both scripts output a file with the combined gene expression counts /example_files/RNA_seq_TPM_files/counts_tpm.txt /example_files/Quant_seq_TPM_files/counts_tpm.txt. This what the /example_files/Quant_seq_TPM_files/counts_tpm.txt looks like:

| gene_id	          | RNA.seq.sample2 | RNA.seq.sample1 | 
| ----------------- | --------------- | --------------- |
| ENSG00000227232.5 | 2.289462046     | 2.13406117      | 
| ENSG00000238009.6 | 0.028186689     | 0.053893268     | 
| ENSG00000241860.6 | 0.143248089     | 0.080565743     | 
| ENSG00000228463.9 |	0.219430279     | 0.173749017     | 

You should also should download or create a clinical.tsv file with details about the cases that will allow you to filter by tumor type or diagnosis. This file can be downloaded from cBioPortal here https://portal.gdc.cancer.gov/repository?files_offset=22700&files_size=100&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22tsv%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%5B%22Gene%20Expression%20Quantification%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%5D%7D&searchTableTab=cases and you can add in your own information for RNA-seq or Quant-seq. The "case_id" column will contain the file name before any extension and you must have the "primary_diagnosis" column as well as that is how samples are grouped, example is here /example_files/clinical.tsv 

# Input gene sets

Gene sets can be imported from the Molecular Signatures Database (MSigDB) or your own gene sets. To use the Molecular Signatures data you must choose the organism and collection which are listed here https://www.gsea-msigdb.org/gsea/msigdb/ an example woudld be "Homo sapiens" "C2". 

`--msigdbr "Homo sapiens" "C2"`

You can also import your own gene sets generated using ChIP_Explorer.sh (https://github.com/chill3456/ChIP_Explorer) which places the gene names in the 3rd to last column example is here /example_files/peaksSet5_all_nearest_gene.bed

`chr1	877011	877838	828	*	0.10530612244897959	0.016902808849117573	chr1	861120	861120	NM_152486	SAMD11	+	15891`

These files can be imported as whole folders. 

`--folderNearestGeneBed "/example_files/gene_sets_folder"`

Files can also be specific as specific files. 

`-specificNearestGeneBed "/example_files/peaksSet4_all_nearest_gene.bed"`

# Running the script 

The script is a R script and can be run on the command line. 

`Rscript /scripts/lsGSEA.R --analysisMethod ssgsea --msigdbr "Homo sapiens" "C2" --folderNearestGeneBed "/example_files/gene_sets_folder" --specificNearestGeneBed "/example_files/peaksSet4_all_nearest_gene.bed" '/example_files/peaksSet5_all_nearest_gene.bed' --cBioCounts "/example_files/cBioPortalFiles/combined_gene_counts.tsv" --RNAseqSampleTpmCounts "/example_files/RNA_seq_TPM_files/counts_tpm.txt" --QuantseqSampleTpmCounts "/example_files/quant_seq_TPM_files/counts_tpm.txt" --outputFolder "/outputExamples/" --clinicalDataSet "/example_files/clinical.tsv"  --setNameToFilter "peaksSet4_all_nearest_gene.bed" "peaksSet1_all_nearest_gene.bed" --primaryDiagnosisFilter "cell_line1" "cell_line2"`

# Arguments

--analysisMethod specifies the type of analysis to run using the GSVA package in R. You can run ssGSEA, GSVA, zscore, or plage. Different methods are described here in the GSVA paper https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7 and the package reference manual https://bioconductor.org/packages/release/bioc/manuals/GSVA/man/GSVA.pdf. ssGSEA is the method that has worked best for data sets such as the example ones. Required argument.

`--analysisMethod ssgsea`

`--analysisMethod GSVA`

`--analysisMethod zscore`

`--analysisMethod plage`

--msigdbr specifies the species and collection if you want to add them from the The Molecular Signatures Database (MSigDB) https://www.gsea-msigdb.org/gsea/msigdb/index.jsp if you want to add them to the gene signatures to test. Optional but must use one of --msigdbr, --analysisMethod, or --specificNearestGeneBed.

`--analysisMethod ssgsea --msigdbr "Homo sapiens" "C2"`

--folderNearestGeneBed specifies a folders containing your own collections of .bed files if you want to add them to the gene signatures to test. The files should contain the gene set in the 3rd to the last column and can be generated using ChIP_Explorer.sh (https://github.com/chill3456/ChIP_Explorer). Optional but must use one of --msigdbr, --analysisMethod, or --specificNearestGeneBed.

`--folderNearestGeneBed "/example_files/gene_sets_folder"`

--specificNearestGeneBed specifies individual .bed files if you want to add them to the gene signatures to test. The files should contain the gene set in the 3rd to the last column and can be generated using ChIP_Explorer.sh (https://github.com/chill3456/ChIP_Explorer). Optional but must use one of --msigdbr, --analysisMethod, or --specificNearestGeneBed.

`--specificNearestGeneBed "/example_files/peaksSet4_all_nearest_gene.bed" '/example_files/peaksSet5_all_nearest_gene.bed'`

--cBioCounts specifies a file containing TPM gene expression counts from the cBioPortal if you want to add them to the samples to run the analysis on. These counts should have been generated using /scripts/expression_matrix_creation_from_cbioportal.py thus have columns gene_id	gene_name	gene_type followed by samples for the other column names. Only the gene_name and sample names columns are used. Optional but use one of --cBioCounts, --analysisMethod, or --specificNearestGeneBed.

`--cBioCounts "/example_files/cBioPortalFiles/combined_gene_counts.tsv"`

--RNAseqSampleTpmCounts specifies a file containing TPM gene expression counts from the RNA-seq if you want to add them to the samples to run the analysis on. These counts should have been generated using /scripts/getting_TPM_counts_from_RNA_seq.R thus have  gene_id column followed by samples for the other column names. Only the gene_id and sample names columns are used. Optional but use one of --cBioCounts, --RNAseqSampleTpmCounts, or --QuantseqSampleTpmCounts.

`--RNAseqSampleTpmCounts "/example_files/RNA_seq_TPM_files/counts_tpm.txt"`

--RNAseqSampleTpmCounts specifies a file containing TPM gene expression counts from the RNA-seq if you want to add them to the samples to run the analysis on. These counts should have been generated using /scripts/getting_TPM_counts_from_RNA_seq.R thus have gene_id column followed by samples for the other column names. Optional but use one of --cBioCounts, --RNAseqSampleTpmCounts, or --QuantseqSampleTpmCounts.

`--RNAseqSampleTpmCounts "/example_files/RNA_seq_TPM_files/counts_tpm.txt"`


--QuantseqSampleTpmCounts specifies a file containing TPM gene expression counts from the RNA-seq if you want to add them to the samples to run the analysis on. These counts should have been generated using /scripts/getting_TPM_counts_from_Quant_seq.R thus have gene_id column followed by samples for the other column names. Optional but use one of --cBioCounts, --RNAseqSampleTpmCounts, or --QuantseqSampleTpmCounts.

`--QuantseqSampleTpmCounts "/example_files/quant_seq_TPM_files/counts_tpm.txt"`

--outputFolder specifies the folder to output the results into. 

`--outputFolder "/outputExamples/"`

--clinicalDataSet specifies the file with details about the cases that will allow you to filter by tumor type or diagnosis. This file can be downloaded from cBioPortal here https://portal.gdc.cancer.gov/repository?files_offset=22700&files_size=100&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22tsv%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%5B%22Gene%20Expression%20Quantification%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%5D%7D&searchTableTab=cases and you can add in your own information for RNA-seq or Quant-seq. The "case_id" column will contain the file name before any extension and you must have the "primary_diagnosis" column as well as that is how samples are grouped, example is here /example_files/clinical.tsv. 

`--clinicalDataSet "/example_files/clinical.tsv"`

--setNameToFilter specifies the gene sets to filter by and plot. Will return a graph showing the average ssGSEA, GSVA, zscore, or plage score and error bars with standard devation for each primary diagnosis.  Will return a .txt file that filtered the file by the given set_name and  a graph showing the    primary_diagnosis's with the 10 lowest and highest average ssGSEA, GSVA, zscore, or plage score and error bars with standard devation.

`--setNameToFilter "peaksSet4_all_nearest_gene.bed" "peaksSet1_all_nearest_gene.bed"`

--primaryDiagnosisFilter specifies the primary_diagnosis to filter by and plot. Will return a .txt file that filtered the file by the given primary_diagnosis and a graph showing the set_name's with the 10 lowest and highest average ssGSEA, GSVA, zscore, or plage score and error bars with standard devation.

`--primaryDiagnosisFilter "cell_line1" "cell_line2"`


# Outputs

A file containing the TPMs for each gene of every sample for being run combined. Is the matrix used for the analysis. /outputExamples/all_TPM_combined.txt

| gene_name | 9bc8ff9e-2f36-5094-ac12-b1bf81b53c15 | 396932ab-c4b3-4526-9201-9222f5f062a5 |	53891f08-a221-4ca8-a502-cf990bdb6020 |	Quant.seq.sample2 |	Quant.seq.sample1 |	RNA.seq.sample2  |	RNA.seq.sample1 |
| --------- | ------------------------------------ | ------------------------------------ | ------------------------------------ | ------------------ | ----------------- | ---------------- | ---------------- |
| A1BG	    | 1.2019                               |	0.1288                              |	5.5751                               |	0.852431894       |	1.086905949       |	0.99503179       |	0.80888827      |
| A1BG-AS1  | 1.4504                               |	1.6928                              |	4.8728                               |	3.488177712       |	2.920839289       |	3.965620003      |	3.813426257     |
| A1CF      | 0                                    |	0.0098                              |	0.0259                               |	0.036894315       |	0.067912375       |	0.002092677      |	0.014071606     |

A file containing the ssGSEA, GSVA, zscore, or plage score for every sample (column) using each gene set (rows) run. The direct output of the analysis /outputExamples/ssgsea_results.txt

|                                | 9bc8ff9e-2f36-5094-ac12-b1bf81b53c15 | 396932ab-c4b3-4526-9201-9222f5f062a5 |	53891f08-a221-4ca8-a502-cf990bdb6020 |	Quant.seq.sample2 |	Quant.seq.sample1 |	RNA.seq.sample2  |	RNA.seq.sample1 |	
| ------------------------------ | ------------------------------------ | ------------------------------------ | ------------------------------------- | ------------------ | ----------------- | ---------------- | ---------------- |
| peaksSet4_all_nearest_gene.bed | 0.073780275                          | 0.066971117                          | -0.031293763                          | 0.028811074        | 0.029973797	      | 0.016757567	     | 0.015351233      |
| peaksSet5_all_nearest_gene.bed | 0.212774822                          | 0.15721897                           | 0.145429982                           | 0.102182261        | 0.10405379 	      |	0.216181993 	   | 0.208751765      |
| peaksSet1_all_nearest_gene.bed | 0.137538644                          | 0.147281865                          | 0.146113581                           | 0.138842354        | 0.137164136 	    |	0.154321421 	   | 0.153168321      |


A file containing the ssGSEA, GSVA, zscore, or plage score (gsva_score column) for every sample (case_id column) using each gene set (set_name column) run merged with the --clinicalDataSet to include the project_id, primary_diagnosis, site_of_resection_or_biopsy, and tissue_or_organ_of_origin columns, and a column saying if it was a postive of negative scrore (up_down column). /outputExamples/ssgsea_merged_results.txt 

| case_id                              | set_name                       |	gsva_score      | project_id |	primary_diagnosis           | site_of_resection_or_biopsy | tissue_or_organ_of_origin | up_down     |
| ------------------------------------ | ------------------------------ | --------------- | ---------- | ---------------------------- | --------------------------- | ------------------------- | ----------- |
| 396932ab-c4b3-4526-9201-9222f5f062a5 | peaksSet4_all_nearest_gene.bed | 0.066971117     | TCGA-LUSC  | Squamous cell carcinoma, NOS | Upper lobe, lung            |	Upper lobe, lung          | Upregulated | 
| 396932ab-c4b3-4526-9201-9222f5f062a5 | peaksSet5_all_nearest_gene.bed | 0.15721897      | TCGA-LUSC	 | Squamous cell carcinoma, NOS | Upper lobe, lung            | Upper lobe, lung          | Upregulated | 
| 396932ab-c4b3-4526-9201-9222f5f062a5 | peaksSet1_all_nearest_gene.bed | 0.147281865     | TCGA-LUSC	 | Squamous cell carcinoma, NOS | Upper lobe, lung            | Upper lobe, lung          | Upregulated | 

A file containing the ssGSEA, GSVA, zscore, or plage score (gsva_score column) averaged for every primary_diagnosis (primary_diagnosis column) by each gene set (set_name column) also returning standard devation (sd column) and and a column saying if it was a postive of negative scrore (up_down column).  /outputExamples/ssgsea_merged_summarized_by_set_name_primary_diagnosis_results.txt

| set_name	                     | primary_diagnosis | gsva_score  | sd          | up_down | 
| ------------------------------ | ----------------- | ----------- | ----------- | ------- |
| peaksSet4_all_nearest_gene.bed | cell_line1        | 0.029392436 | 0.00082217  | red     | 
| peaksSet4_all_nearest_gene.bed | cell_line2        | 0.0160544	 | 0.000994428 | red     | 
| peaksSet5_all_nearest_gene.bed | cell_line1        | 0.103118025 | 0.001323371 | red     | 

If using --setNameToFilter argument will return a .txt file that filtered the /outputExamples/ssgsea_merged_summarized_by_set_name_primary_diagnosis_results.txt file by the given set_name and a graph showing the    primary_diagnosis's with the 10 lowest and highest average ssGSEA, GSVA, zscore, or plage score and error bars with standard devation.  /outputExamples/peaksSet4_all_nearest_gene.bed__filtered_merged_summarized_by_set_name_primary_diagnosis_results.txt and /outputExamples/peaksSet4_all_nearest_gene.bed__filtered_merged_summarized_by_set_name_primary_diagnosis_results.pdf

| primary_diagnosis                       | gsva_score  | sd          | up_down | 
| --------------------------------------- | ----------- | ----------- | ------- | 
| Neuroblastoma, NOS                      | 0.073780275	| NA          | red     | 
| Squamous cell carcinoma, NOS            | 0.066971117 | NA	        | red     | 
| cell_line1                              | 0.029392436	| 0.00082217  | red     | 
| cell_line2                              | 0.0160544   | 0.000994428	| red     | 
| Precursor B-cell lymphoblastic leukemia	| -0.031293763| NA          | blue    | 

<img src='https://github.com/chill3456/lsGSEA/blob/master/assets/peaksSet4_all_nearest_gene.bed__filtered_merged_summarized_by_set_name_primary_diagnosis_results.png' width = '50%' alt="graph showing peaksSet4_all_nearest_gene.bed__filtered_merged_summarized_by_set_name_primary_diagnosis_results.png">

If using --primaryDiagnosisFilter argument will return a .txt file that filtered the /outputExamples/ssgsea_merged_summarized_by_set_name_primary_diagnosis_results.txt file by the given set_name and a graph showing the    set_name's with the 10 lowest and highest average ssGSEA, GSVA, zscore, or plage score and error bars with standard devation.  /outputExamples/cell_line1__filtered_merged_summarized_by_set_name_primary_diagnosis_results.txt and /outputExamples/cell_line1__filtered_merged_summarized_by_set_name_primary_diagnosis_results.pdf

| set_name                       | gsva_score  | sd          | up_down | 
| ------------------------------ | ----------- | ----------- | ------- | 
| peaksSet4_all_nearest_gene.bed | 0.029392436 | 0.00082217  | red     | 
| peaksSet5_all_nearest_gene.bed | 0.103118025 | 0.001323371 | red     | 
| peaksSet1_all_nearest_gene.bed | 0.138003245 | 0.001186679 | red     | 
| peaksSet2_all_nearest_gene.bed | 0.25546955  | 0.008055932 | red     | 
| peaksSet3_all_nearest_gene.bed | 0.143625617 | 0.0022123   | red     | 
| ABBUD_LIF_SIGNALING_1_DN       | 0.030325599 | 0.003944133 | red     | 

<img src='https://github.com/chill3456/lsGSEA/blob/master/assets/cell_line1__filtered_merged_summarized_by_set_name_primary_diagnosis_results.png' width = '50%' alt="graph showing cell_line1__filtered_merged_summarized_by_set_name_primary_diagnosis_results.png">
