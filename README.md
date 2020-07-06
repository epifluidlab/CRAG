# CRAG
**C**ell f**R**ee dn**A** fra**G**mentation (CRAG): *De novo* characterization of cell-free DNA fragmentation hotspots.

## Table of Contents
1. [Quick start](#Quick-start)
2. [Installation](#Installation)
3. [Usage](#Usage)
4. [Citation](#Citation)
5. [License](#License)
6. [Contact](#Contact)

## Quick start
1. Install
	```
	git clone --recursive https://github.com/epifluidlab/CRAG.git
	cd CRAG
  	pip install pysam
	```

2. Run the example test file at your bash command line to call hotspot
	```
	cd CRAG/hotspots_calling
	python bam_read.py -in test.bam -out test_dir
	matlab -nodisplay -r 'CRAG test_dir 1; exit;' 	
	```
	
	This should produce the following files (inside test_dir/result_n/):
	* the hotspots (peak_all.mat, peaks.bed), 
	* fragmentation pattern around hotspots (IFS_plot.fig, IFS_plot.pdf)
	* enrichment patterns (TSS_plot.fig, TSS_plot.pdf, CTCF_plot.fig, CTCF_plot.pdf, Enrichment_plot.fig,Enrichment_plot.pdf) 

## Installation
#### Prerequisites
* Linux, Mac OSX, and Windows
* Matlab 2019b
* python 2.7 
* pysam 0.12.0.1 or above
* samtools 1.9

#### required files
1. indexed Bam file (paired-end whole-genome sequencing, recommend to have at least 400 million fragments in autosomes after the samtools filtering step. If you want to call hotspots for several chrommsomes (not the whole autosome), you can provide the bam file only with the corresponding chromosomes.)
2. Dark regions (provided in resource folder for hg19, or could be downloaded from UCSC table browser)
3. Mappability files (provided in resource folder for hg19, or could be generated by [GEM library tools](https://wiki.bits.vib.be/index.php/Create_a_mappability_track)

## Usage

### Sample pre-processing and read
```
samtools view -bh -f 3 -F 3852 -q 30 input.bam > output.filtered.bam

python bam_read.py -in output.filtered.bam -out result_dir
```

### Hotspot calling
Two modes for the hotspot calling: 1 - call hotspots using IFS. 2 - call hotspots using GC bias corrected IFS.
There are two choices to run the tool:
- Run it directly at matlab (all the following code examples will follow this style)
```
CRAG('result_dir',1)
```
- Run it at linux bash with matlab installed
```
matlab -nodisplay -r 'CRAG result_dir 1; exit;'
```
#### Options for hotspot calling

* global p value cut-off: argument: 'global_p', a positive number between 0 to 1. default: 1e-5. 
* local p value cut-off: argument: 'local_p', a positive number between 0 to 1.default: 1e-5
* fdr cut-off:argument: 'fdr', a positive number between 0 to 1. default: 0.01
* hotspot distance for merging: argument: 'distance', a positive integer. default: 200
* whether or not do enrichment for the hotspots:  argument: 'enrichment', 0 or 1. default: 1. 

Example code with option specified:
`CRAG('result_dir',1,'local_p',0.001,'distance',300,'enrichment',0)`

This command will call hotspots in 'result_dir', using IFS without GC-bias correction, with global p value cut-off = 0.00001, local p value cut-off = 0.001 and fdr cut-off = 0.01. After the significant regions were detected, the regions nearby (with distance less than 300bp) were merged. And the pipeline will produce four files (two hotspot files: hotspots.bed and hotspots.mat and two files of the fragmentation pattern around hotspots (IFS_plot.fig, IFS_plot.pdf))in result_dir/result/n.

### Unsupervised clustering, PCA, and t-SNE
You will need to run the script 'bam_read.py' to read the fragment for all the samples and write them to txt files.

	`cd CRAG/unsupervised_analysis`
	
1. Obtain IFS score in each sample. If using low-coverage WGS, please group samples based on their category. e.g.:

	`IFS_write('PNAS/C309')`
	
2. Call hotspots for all the samples of a specific categoty:
	* Function:
		`Hotspot_call_multi_sample.m`
	* Parameters:
		* file_list: an file name of excel (.xlsx or xls) file, which contain all the sample names, or a cell variable (N /* 1) contain all the sample names. For example, 'HCC.xlsx', which contain one column, in each row of the column is the name a sample of HCC.
		* input_path:the path of all the samples, such as 'PNAS'
		* out_name:the outfile file name you want, such as 'HCC'
		* peak_type: 1-call hotspots without GC bias correction; 2- call hotspots based on IFS after GC bias correction
	* Optional parameters (the same with hotspot calling in single sample):
		* global_p: global p-value cut off
		* local_p: p-value cut-off for local test
		* fdr: cut-off
		* distance: Distance cut-off to merge the significant regions nearby.
		* enrichment:whether or not do enrichment for the hotspot.
	* Example code to call hotspot (IFS, no GC bias corrected) in all HCC samples under PNAS/ directory by using the sample information at HCC.xlsx:
	`Hotspot_call_multi_sample('HCC.xlsx','PNAS','HCC',1,'enrichment',0)`
	
3. Get IFS matrix for each hotspot in each sample: get the hotspots' IFS of each sample and merge all the samples of a specific category as a matrix and output the matrix:
	* Function: 
		`IFS_matrix_obtain.m`
	* Parameters:
		* file_list:the name of an excel file (include the Suffix name), or a variable include all the the files.
		* input_path: The path of all the files
		* out_name: The output file name (matrix file name)
		* peak_type： Peak_type==1: IFS. Peak_type==2: GC bias corrected IFS.
		* peak_file: The sample name for the hotspots to obtain the IFS score.
	* For example, when we do upsupervised analysis of HCC samples and healthy samples, For each HCC samples (healthy samples), we want to get the IFS score of both the HCC hotspots and healthy hotspots. If we have the hotspots of HCC samples in HCC/result_n/ and the hotspots of healthy samples in healthy/result_n/, we could run:
	`IFS_matrix_obtain('HCC.xlsx','PNAS','HCC/result_n/norm_IFS.mat',1,'HCC','healthy')`
	
4. Unsupervised hierarchical clustering: do hierarchical clustering for several types of samples based on the IFS matrix (selected features by variance).
	* Function: 
		`cluster_analysis_var.m`
	* Parameters:
		* distance_type:
			* 1: clustergram using spearman Correlation coefficient as similarity (distance evalulation method) and 'weighted' as linkage algorithm;
			* 2: clustergram using euclidean as similarity (distance evalulation method) and 'ward' as linkage algorithm
		* feature_num: The number of the hotspots were used for clustering (most variabel hotpots).
		* matrix_file：The .mat file for the IFS matrix of this type 
		* color:after clustering, the color to denote this type

	* For example, if we want to use top 10k hotspots in HCC and healthy samples to do hierarchical clustering (Spearman Correlation coefficient as similarity (distance evalulation method) and 'weighted' as linkage algorithm). "red" to color HCC and "green" to color healthy samples. The command should be:
		`cluster_analysis_var(1,100000,'HCC/result_n/norm_IFS.mat','r','healthy/result_n/norm_IFS.mat','g')`
	
5. PCA analysis:
	* Function: 
		`PCA_analysis.m`
	* Parameters:
		* color:after clustering, the color to denote this type
	* Example code:
		`PCA_analysis('HCC/result_n/norm_IFS.mat','r','healthy/result_n/norm_IFS.mat','g')`
	
6. One-way anvoa analysis to select the most variable features:
	* Function: 
		`PCA_analysis.m`
	* Parameters:
		* p_cut: The p-value cut off in one-way anvoa to select the features
		* out_name: The name of the output data. 
		* matrix_file1：The .mat file for the IFS matrix of the first type
		* color1: The color for the first category (for future analysis).
		* matrix_file2：The .mat file for the IFS matrix of the second type
		* color2: The color for the second category.
		* matrix_fileN
		* colorN
	* For example, if we want to visualize the breast cancer samples (IFS matrtrix in Breast/norm_IFS.mat), lung cancer samples (IFS matrtrix in Lung/norm_IFS.mat), ovarian cancer (IFS matrtrix in Ovarian/norm_IFS.mat) and healthy samples (IFS matrtrix in healthy/norm_IFS.mat), the command could be as following:
		`ANVOA_calculate (0.01,'ANVOA_IFS_matrix.mat','Breast/norm_IFS.mat','m','Lung/norm_IFS.mat','k','Lung/norm_IFS.mat','b','healthy/norm_IFS.mat','r')`
	After running this code, we could get a .mat file (ANVOA_IFS_matrix.mat) in the current workplace, in the file, it contains four variables:
		* data: M * N, N is the sample number of the four types and M is the hotspots that have p-value in one-way anvoa no less than 0.01.
		* label:N * 1, the label of all the samples (Breast cancer sample as 1, Lung cancer sample as 2...)
		* p: M * 1,p-value of the selected hotspots
		* col:{'m';'k','b','r'}, the color of the four types, which could be used in TSNE.
		
7. run t-SNE:
	* Function: 
		`TSNE.m`
	* Parameters:
		* matrix file: The file should including at least three variable:
			* data: N * M matrix,N - hotspots number,M - sample num;
			* label:M*1: The category of the M samples (1,2,3,etc)
			* p：N*1：p-value of the hotspots in anvoa test
		* col:X * 1,cell of string, describes the color of the X types, if not exist, the color will be randomly assigned.
		* distance_type: 
			* 1: the distance among the samples will be evalulate by spearman Correlation coefficient
			* 2: the distance among the samples will be evalulate by euclidean distance
			* neighbor_num: the number of neighbors (similar as plexity in python/R)
			* p_cut, the cut-off to choose the features by one-way ANOVA test (i.e. 0.01).
	* Example code:
		`TSNE(1,30,0.01,'ANVOA_IFS_matrix.mat')`
	Spearman Correlation coefficient will be used to calculate the distance among the samples, 30 neighbors will be used in TSNE analysis, and all the the hotspots with anvoa p-value <= 0.01 will be used in TSNE.
	* Output:
		* two figures (TSNE_plot.pdf and TSNE_plot.fig) 
		* a .mat file （tsne_score.mat） contains two variables:
			* class_label：The labels of the samples
			* score: The scores of the first two components of each samples.
	
	
8. Unsupervised clustering after selecting features by one-way ANOVA test:
	`cluster_nature_GC()`


### Cancer vs. healthy classification
TBD

#### Options for cancer vs. healthy classification
```
- k fold
- random seed
- GC bias correction
```

* Recommanded setting:
`code`

### Tissues-of-origin prediction
TBD

#### Options for tissues-of-origin prediction
```
- k fold
- random seed
- GC bias correction
```

* Recommanded setting:
`code`

#### Raw source code in manuscript
The raw source code and their readme files for the publication are inside manuscript/

## Citation
```
Zhou X & Liu Y. "De novo characterization of cell-free DNA fragmentation hotspots boosts the power for early detection and localization of multi-cancer". In preparation.
```
## License
* For academic research, please refer to MIT license.
* For commerical usage, please contact the authors.

## Contact
* Xionghui Zhou <xionghui.zhou@cchmc.org>
* Yaping Liu <lyping1986@gmail.com>















	

