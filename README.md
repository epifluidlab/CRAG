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
2. Download the test bam file and required GC and mappability files
	```
	cd Basic_info
	wget -c https://zenodo.org/record/3928546/files/GC.zip
	wget -c https://zenodo.org/record/3928546/files/mappability.zip
	wget -c https://zenodo.org/record/3928546/files/chrm_state.zip
	unzip GC.zip
	unzip mappability.zip
	unzip chrm_state.zip
	cd ../hotspots_calling
	wget -c https://zenodo.org/record/3928546/files/BH01.chr22.bam
	wget -c https://zenodo.org/record/3928546/files/BH01.chr22.bam.bai
	```
3. Run the example test file at your bash command line to call hotspot. (**link the Basic_info/ into your current working directory**)
	```
	ln -s ../Basic_info
	python bam_read.py -in BH01.chr22.bam -out test_dir
	matlab -nodisplay -r 'CRAG test_dir 1; exit;' 	
	```
**You need at least 10Gb memory to finish the test example**. At our server, it costs about 10 mins at CentOS 7 with 10Gb memory and one CPU core: Intel(R) Xeon(R) CPU E5-2695 v3 @ 2.30GHz	

This should produce the following files (inside test_dir/result_n/):
	* the hotspots (peak_all.mat, peaks.bed), 
	* fragmentation pattern around hotspots (IFS_plot.fig, IFS_plot.pdf)
	* enrichment patterns (TSS_plot.fig, TSS_plot.pdf, CTCF_plot.fig, CTCF_plot.pdf, Enrichment_plot.fig, Enrichment_plot.pdf)  

## Installation
#### Prerequisites
* Linux, Mac OSX, and Windows (with at least 10Gb memory for each CPU core)
* Matlab 2019b
* python 2.7 
* pysam 0.12.0.1 or above
* samtools 1.9


#### required files
* Indexed Bam file (paired-end whole-genome sequencing, recommend to have at least 400 million fragments in autosomes after the samtools filtering step. If you want to call hotspots for several chrommsomes (not the whole autosome), you can provide the bam file only with the corresponding chromosomes.)
* Basic_info directory **Always link the Basic_info/ under your current working directory** (example is showed in Quick Start part)
* GC content files (provided in [zenodo.org](https://zenodo.org/record/3928546/files/GC.zip) for hg19/GRch37, download it under Basic_info directory and unzip it)
* Mappability files (provided in [zenodo.org](https://zenodo.org/record/3928546/files/mappability.zip) for hg19/GRch37, download it under Basic_info directory and unzip it)
* ChromHMM state files (provided in [zenodo.org](https://zenodo.org/record/3928546/files/chrm_state.zip) for hg19/GRch37, download it under Basic_info directory and unzip it)

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
	```
	cd CRAG/unsupervised_analysis
	```
	
1. Obtain IFS score in each sample. If using low-coverage WGS, please group samples based on their category. e.g.:
	```
	IFS_write('PNAS/C309')
	```
	
2. Call hotspots for all the samples of a specific categoty:
	* Function:
		```
		Hotspot_call_multi_sample.m
		```
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
		```
		Hotspot_call_multi_sample('HCC.xlsx','PNAS','HCC',1,'enrichment',0)
		```
	
3. Get IFS matrix for each hotspot in each sample: get the hotspots' IFS of each sample and merge all the samples of a specific category as a matrix and output the matrix:
	* Function: 
		```
		IFS_matrix_obtain.m
		```
	* Parameters:
		* file_list:the name of an excel file (include the Suffix name), or a variable include all the the files.
		* input_path: The path of all the files
		* out_name: The output file name (matrix file name)
		* peak_type： Peak_type==1: IFS. Peak_type==2: GC bias corrected IFS.
		* peak_file: The sample name for the hotspots to obtain the IFS score.
	* For example, when we do upsupervised analysis of HCC samples and healthy samples, For each HCC samples (healthy samples), we want to get the IFS score of both the HCC hotspots and healthy hotspots. If we have the hotspots of HCC samples in HCC/result_n/ and the hotspots of healthy samples in healthy/result_n/, we could run:
		```
		IFS_matrix_obtain('HCC.xlsx','PNAS','HCC/result_n/norm_IFS.mat',1,'HCC','healthy')
		```
	
4. Unsupervised hierarchical clustering: do hierarchical clustering for several types of samples based on the IFS matrix (selected features by variance).
	* Function: 
		```
		cluster_analysis_var.m
		```
	* Parameters:
		* distance_type:
			* 1: clustergram using spearman Correlation coefficient as similarity (distance evalulation method) and 'weighted' as linkage algorithm;
			* 2: clustergram using euclidean as similarity (distance evalulation method) and 'ward' as linkage algorithm
		* feature_num: The number of the hotspots were used for clustering (most variabel hotpots).
		* matrix_file：The .mat file for the IFS matrix of this type 
		* color:after clustering, the color to denote this type

	* For example, if we want to use top 10k hotspots in HCC and healthy samples to do hierarchical clustering (Spearman Correlation coefficient as similarity (distance evalulation method) and 'weighted' as linkage algorithm). "red" to color HCC and "green" to color healthy samples. The command should be:
		```cluster_analysis_var(1,100000,'HCC/result_n/norm_IFS.mat','r','healthy/result_n/norm_IFS.mat','g')
		```
	
5. PCA analysis:
	* Function: 
		```
		PCA_analysis.m
		```
	* Parameters:
		* color:after clustering, the color to denote this type
	* Example code:
		```
		PCA_analysis('HCC/result_n/norm_IFS.mat','r','healthy/result_n/norm_IFS.mat','g')
		```
	
6. One-way anvoa analysis to select the most variable features:
	* Function: 
		```
		PCA_analysis.m
		```
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
		```
		ANVOA_calculate (0.01,'ANVOA_IFS_matrix.mat','Breast/norm_IFS.mat','m','Lung/norm_IFS.mat','k','Lung/norm_IFS.mat','b','healthy/norm_IFS.mat','r')
		```
	* Output:
		* a .mat file (ANVOA_IFS_matrix.mat) in the current workplace, in the file, it contains four variables:
			* data: M * N, N is the sample number of the four types and M is the hotspots that have p-value in one-way anvoa no less than 0.01.
			* label:N * 1, the label of all the samples (Breast cancer sample as 1, Lung cancer sample as 2...)
			* p: M * 1,p-value of the selected hotspots
			* col:{'m';'k','b','r'}, the color of the four types, which could be used in TSNE.
		
7. run t-SNE:
	* Function: 
		```
		TSNE.m
		```
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
		```
		TSNE(1,30,0.01,'ANVOA_IFS_matrix.mat')
		```
	Spearman Correlation coefficient will be used to calculate the distance among the samples, 30 neighbors will be used in TSNE analysis, and all the the hotspots with anvoa p-value <= 0.01 will be used in TSNE.
	* Output:
		* two figures (TSNE_plot.pdf and TSNE_plot.fig) 
		* a .mat file （tsne_score.mat） contains two variables:
			* class_label：The labels of the samples
			* score: The scores of the first two components of each samples.
	
	
8. Unsupervised clustering after selecting features by one-way ANOVA test:
	* Function: 
		```
		cluster_analysis_ANVOA.m
		```
	* Parameters:
		* matrix file: The file should including at least three variable:
			* data:N * M matrix,N - hotspots number,M - sample num;
			* label:M * 1: The category of the M samples (1,2,3,etc)
			* p：N * 1：p-value of the hotspots in anvoa test
			* col:X * 1,cell of string, describes the color of the X types, if not exist, the color will be randomly assigned.
			* distance_type:
				* 1: clustergram using spearman Correlation coefficient as similarity (distance evalulation method) and 'weighted' as linkage algorithm;
				* 2: clustergram using euclidean as similarity (distance evalulation method) and 'ward' as linkage algorithm.
			* p_cut:the cut-off to choose the features (i.e. 0.01).
			* feature_num:The number of the type-specific hotspots in each type.
	* Example code:
		```
		cluster_analysis_ANVOA(1,0.01,5000,'ANVOA_IFS_matrix.mat')
		```
	This command will run on the anvoa filtered matrix, p-value <= 0.01 will be used, and for each category of the data,
5000 hotspots will be selected by the z-score difference, based on these hotspots, clustergram, with spearman Correlation coefficient as similarity (distance evalulation method) and 'weighted' as linkage algorithm, will be used to do hierarchical clustering.
	* Output:
		* two figures (clustering_result.fig and clustering_result.pdf)
		* one .mat file contain the matrix that was used to do clustering analysis and the category information of these samples (cluster_anvoa_data.mat)

### Cancer vs. healthy classification and tissues-of-origin prediction
You should run the script 'bam_read.py' to read the fragment for all the samples and write it to txt files.We recommend put all the samples from the same data set into the same folder. For example, when we analysis the Sun et. al. data set, we should put the fragment information of all the samples in the folder PNAS
```
cd CRAG/classification
```
1. Obtain IFS score in each sample. If using low-coverage WGS, please group samples based on their category. e.g.:
	```
	IFS_write('PNAS/C309')
	```
2. Call hotspots for each category only in training fold in k-fold cross validation. If do 10-fold validation, the program should be executed 10 times, change ith from 1 to 10:
	* Function:
		```
		k_fold_call.m
		```
	* Parameters:
		* fold_info: an file name of excel (.xlsx or xls) file, which contain all the sample names (string),category (string) and fold id (interger,i.e.1,2,3...), or a cell variable (N * 3) contain all this information. For example, 'HCC_healthy_fold.xlsx', which contain three column, in each row,the first column is its sample name, the second column is its category ('HCC' or 'healthy') and the last column is the fold id it begongs to (i.e. 5).
		* input_path: the path that the fragment information of these samples located in (i.e.PNAS).
		* out_name_name: The folder the all the output should be saved (i.e. HCC_healthy).
		* peak_type: 1-call hotspots without GC bias correction; 2- call hotspots based on IFS after GC bias correction
		* ith: call hotspots for which fold in the cross validation (i.e. 2).
	* Optional parameters (the same with hotspot calling in single sample):
		* global_p: global p-value cut off
		* local_p: p-value cut-off for local test
		* fdr: cut-off
		* distance: Distance cut-off to merge the significant regions nearby.
		* enrichment:whether or not do enrichment for the hotspots.		
	* Example code to call hotspot (IFS, no GC bias corrected) in all HCC samples under PNAS/ directory by using the sample information at HCC.xlsx:
		```
		k_fold_call('HCC_healthy_fold.xlsx','PNAS','HCC_healthy',1,2,'enrichment',0)
		```
	* Output:
		After run the function, in the folder of 'HCC_healthy', the program will produce a file 'HCC_healthy/2/HCC/, and in HCC_healthy/2/HCC/result_n/, four files would exist (peak_all.mat,peaks.bed) and two files of the fragmentation pattern around hotspots (IFS_plot.fig, IFS_plot.pdf)). Also, the program will produce a file 'HCC_healthy/5/healthy/, and in HCC_healthy/5/healthy/result_n/, four files would exist (peak_all.mat,peaks.bed) and two files of the fragmentation pattern around hotspots (IFS_plot.fig, IFS_plot.pdf))

3. Get IFS matrix in each fold: get the hotspots' IFS of each sample and merge all the samples of a specific category in the training set as a matrix and output the matrix, also the program will output the matrix for the testing samples. If do 10-fold validation, the program should be executed 10 times, change ith from 1 to 10:
	* Function:
		```
		k_fold_IFS_matrix.m
		```
	* Parameters:
		* fold_info: an file name of excel (.xlsx or xls) file, which contain all the sample names (string),category (string) and fold id (interger,i.e.1,2,3...), or a cell variable (N * 3) contain all this information. For example, 'HCC_healthy_fold.xlsx', which contain three column, in each row,the first column is its sample name, the second column is its category ('HCC' or 'healthy') and the last column is the fold id it begongs to (i.e. 5).
		* input_path: the path that the fragment information of these samples located in (i.e.PNAS).
		* out_name_name: The folder the all the output should be saved (i.e. HCC_healthy).
		* peak_type: 1-get IFS of the hotspots without GC bias correction; 2- get IFS of hotspots based on IFS after GC bias correctionand correct the IFS.
		* ith: ontain the IFS matrix for which fold in the cross validation (i.e. 2).
		* pos_class:a string shows the positive class in the classification, i.e. 'HCC', if this is a fold-validation for multi-classification,spefiific it as 'None'.
	* Example code:
		```
		k_fold_IFS_matrix('HCC_healthy_fold.xlsx','PNAS','HCC_healthy',1,2,'HCC')
		```
	* Output:
		In the folder of 'HCC_healthy', the program will produce files in 'HCC_healthy/2/result_n/. Three files would exist (train_norm_data1.mat,train_norm_data2.mat,test_norm_data.mat)
			* In train_norm_data1.mat, there is a M * N matrix, N is the number of HCC samples in the second training set. M is the total hotspot number of HCC samples and healthy samples in this fold.
			* In train_norm_data2.mat, there is a M * N2 matrix, N2 is the number of Healthy samples in the second training set. M is the total hotspot number of HCC samples and healthy samples in this fold.
			* In test_norm_data.mat, there is a M * N3 matrix, N3 is the number of samples in the second test data set. M is the total hotspot number of HCC samples and healthy samples in this fold (training set). Also, there is a variable 'test_label' (N3 * 1), is the class labels of all the test samples.

4. For multi-class classification in tissues-of-origin analysis:
For example, 'Breast cancer vs. lung cancer vs. ovarian cancer vs. colon cancer vs. Pancreatic cancer', we only need to change a parameter in the function.
e.g.
	* `k_fold_call('multi_cancer_fold.xlsx','multi_cancer','TOO',2,2,'enrichment',0)`
	* `k_fold_IFS_matrix('multi_cancer_fold.xlsx','multi_cancer','TOO',2,2,'None')`

5. Binary classification (e.g. Cancer vs. healthy) by liner svm model:
	* Function:
		```
		IFS_binaryClass_classification.m
		```
	* Parameters:
		* input_file_path: The file of the matrix file, e.g. 'HCC_healthy'
		* fold_number: How many folds are in the validation, e.g. 10
		* feature_num: The number of features used.
	* Example code:
		```
		IFS_binaryClass_classification('HCC_healthy',10,30000)
		```
	After this command, the program would perform 10-fold validation classification in the folder 'HCC_healthy', using the most stable 30,000 hostpots as features, and liner svm as model. At last, in the floder 'HCC_healthy/', there would be two figures (ROC.fig and ROC.PDF), also there would be a .mat file (classfication_result.mat), which contains all the detail information in the classification.

6. Multi-class classification (e.g. tissues-of-origin):
	* Function:
		```
		IFS_binaryClass_classification.m
		```
	* Parameters:
		* input_file_path: The file of the matrix file, e.g. 'TOO'
		* fold_number: How many folds are in the validation, e.g. 10
		* class_num: Howm many categories are in this classification, e.g. 5
		* feature_num: The number of features used (for the step in decision tree).
	* Example code:
		```
		IFS_multiClass_classification('TOO',10,5,100000)
		```
	After this command, the program would perform 10-fold validation classification in the folder 'TOO', using centroid distance to do top 2 predictions, and then the most stable 100,000 hostpots as features, and decision tree as model were used to decide top 1 prediction. At last, in the floder 'TOO/', there would be a .mat file (classfication_result.mat), which contains all the detail information in the classification:
		* The furst two columns: top 2 predictor predicted by centroid distance
		* The third column: The top 1 predictor predicted by decision tree
		* The forth column: The true class label
		* The fifth column: The fold id

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















	

