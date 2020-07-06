(1)The folder contains all the code for upsupervised analysis, including clustergram, TSNE and PCA analysis.
(2)To run this pipeline, we should run the script 'bam_read.py' to read the fragment for all the samples and write it to txt files.  
For example:python bam_read.py -in C309.filter.bam -out PNAS/C309 could be used to write all the fragment information into the path PNAS/C309. We recommend put all the samples from the same data set into the same folder. For example, when we analysis the Sun et. al. data set, we should put the fragment information of all the samples in the folder PNAS.
(3)After that, IFS_write.m could be used to get the IFS socre for all the chromosomes of the sample. 
For example:IFS_write('PNAS/C309').
(4)Hotspot_call_multi_sample.m could be used to call hotspots for all the samples of a specific categoty, for example, all the HCC samples in Sun et. al. data set.
In the function, the parameters are shown as follow:
file_list: an file name of excel (.xlsx or xls) file, which contain all the sample names, or a cell variable (N*1) contain all the sample names. For example, 'HCC.xlsx', which contain one column, in each row of the column is the name a sample of HCC.
input_path:the path of all the samples, such as 'PNAS'
out_name:the outfile file name you want, such as 'HCC'
peak_type: 1-call hotspots without GC bias correction; 2- call hotspots based on IFS after GC bias correction
optional parameters (the same with hotspot calling in single sample):
global_p: global p-value cut off
local_p: p-value cut-off for local test
fdr: cut-off
distance: Distance cut-off to merge the significant regions nearby.
enrichment:whether or not do enrichment for the hotspots:  argument: '
An example to run the pipeline is: Hotspot_call_multi_sample('HCC.xlsx','PNAS','HCC',1,'enrichment',0);
After run the function, the program will produce a file HCC, and in HCC/result_n/, four files would exist (peak_all.mat,peaks.bed) and two files of the fragmentation pattern around hotspots (IFS_plot.fig, IFS_plot.pdf)) in HCC/result_n.
(5)IFS_matrix_obtain.m could be used to get the hotspots' IFS of each sample and merge all the samples of a specific category as a matrix and output the matrix.
In the function, the parameters are shown as follow:
file_list:the name of an excel file (include the Suffix name), or a variable include all the the files
input_path: The path of all the files
out_name: The output file name (matrix file name)
peak_type： Peak_type==1: IFS. Peak_type==2: GC bias corrected IFS.
peak_file: The sample name for the hotspots to obtain the IFS score.
optional parameters:
Several sample names of several hotspot sets
For example, when we do upsupervised analysis of HCC samples and healthy samples, For each HCC samples (healthy samples), we want to get the IFS score of both the HCC hotspots and healthy hotspots. If we have the hotspots of HCC samples in HCC/result_n/ and the hotspots of healthy samples in healthy/result_n/, we could run:
IFS_matrix_obtain('HCC.xlsx','PNAS','HCC/result_n/norm_IFS.mat',1,'HCC','healthy').
Also, IFS_matrix_obtain('healthy.xlsx','PNAS','healthy/result_n/norm_IFS.mat',1,'HCC','healthy') could be used to get the IFS matrix for healthy samples.
(6)cluster_analysis_var.m could be used to do hierarchical clustering for several types of samples based on the IFS matrix (selected features by variance).
In the function, the parameters are shown as follow:
distance_type =1, clustergram using spearman Correlation coefficient as similarity (distance evalulation method) and 'weighted' as linkage algorithm
distance_type =2:clustergram using euclidean as similarity (distance evalulation method) and 'ward' as linkage algorithm
feature_num: The number of the hotspots were used for clustering (most variabel hotpots).
matrix_file：The .mat file for the IFS matrix of this type 
color:after clustering, the color to denote this type
Options: 2*N (N>=1) variable, the matrix file and its corresponding color
For example, if we want to use top 10k hotspots in HCC and healthy samples to do hierarchical clustering (Spearman Correlation coefficient as similarity (distance evalulation method) and 'weighted' as linkage algorithm).
The command should be: cluster_analysis_var(1,100000,'HCC/result_n/norm_IFS.mat','r','healthy/result_n/norm_IFS.mat','g');
After that, in the current workplace, there would two figures (clustering_result.fig and clustering_result.pdf) and one .mat file contain the matrix that was used to do clustering analysis and the category information of these samples.
(7)For PCA analysis, taking PCA analysis in HCC vs. healthy samples as an example.
PCA_analysis.m could be used.
In the function, the parameters are shown as follow:
color:after clustering, the color to denote this type
Options: 2*N (N>=1) variable, the matrix file and its corresponding color
The command could be like this: PCA_analysis('HCC/result_n/norm_IFS.mat','r','healthy/result_n/norm_IFS.mat','g');
After that, in the current workplace, there would two figures (PCA_plot.fig and PCA_plot.pdf) and one .mat file contain the sccors of all the samples after PCA analysis and the category information of these samples.
(8)For TSNE analysis, the ANVOA_calculate.m should be used to do one-way anvoa analysis for the IFS matrix to select the features.
In the function, the parameters are shown as follow:
p_cut: The p-value cut off in one-way anvoa to select the features
out_name: The name of the output data. 
matrix_file1：The .mat file for the IFS matrix of the first type
color1: The color for the first category (for future analysis).
matrix_file2：The .mat file for the IFS matrix of the second type
color2: The color for the second category.
Options:more the matrix files and the coressponding color.
For example, if we want to visualize the breast cancer samples (IFS matrtrix in Breast/norm_IFS.mat), lung cancer samples (IFS matrtrix in Lung/norm_IFS.mat), ovarian cancer (IFS matrtrix in Ovarian/norm_IFS.mat) and healthy samples (IFS matrtrix in healthy/norm_IFS.mat), the command could be as following:
ANVOA_calculate (0.01,'ANVOA_IFS_matrix.mat','Breast/norm_IFS.mat','m','Lung/norm_IFS.mat','k','Lung/norm_IFS.mat','b','healthy/norm_IFS.mat','r');
After running this code, we could get a .mat file (ANVOA_IFS_matrix.mat) in the current workplace, in the file, it contains four variables:
data: M*N, N is the sample number of the four types and M is the hotspots that have p-value in one-way anvoa no less than 0.01.
label:N*1, the label of all the samples (Breast cancer sample as 1, Lung cancer sample as 2...)
p: M*1,p-value of the selected hotspots
col:{'m';'k','b','r'}, the color of the four types, which could be used in TSNE.
(9)TSNE anlaysis:
Function: TSNE(distance_type,neighbor_num,p_cut,matrix_file)
%%%%%%%%TSNE visualize the samples from the matrix file after anvoa analysis.
matrix file: The file should including at least three variable:
data:N*M matrix,N - hotspots number,M - sample num;
label:M*1: The category of the M samples (1,2,3,etc)
p：N*1：p-value of the hotspots in anvoa test
optional:col:X*1,cell of string, describes the color of the X types, if not exist, the color will be randomly assigned.
distance_type =1, the distance among the samples will be evalulate by spearman Correlation coefficient
distance_type =2: the distance among the samples will be evalulate by euclidean distance
neighbor_num: the number of neighbors
p_cut, the cut-off to choose the features (i.e. 0.01).
e.g. TSNE(1,30,0.01,'ANVOA_IFS_matrix.mat') could do TSNE analysis for the matrix produced above. Spearman Correlation coefficient could be used to calculate the distance among the samples, 30 neighbors will be used in TSNE analysis, and all the the hotspots with anvoa p-value <= 0.01 would be used in TSNE.
The output of the command would be two figures (TSNE_plot.pdf and TSNE_plot.fig) and a .mat file （tsne_score.mat） contains two variables:
class_label：The labels of the samples
score: The scores of the first two components of each samples.
(10)cluster_analysis_ANVOA.m could be used to do hierarchical clustering for several types of samples based on the selected features by ANVOA.
Function: cluster_analysis_ANVOA(distance_type,p_cut,feature_num,matrix_file)
In the function, the parameters are shown as follow:

matrix file: The file should including at least three variable:
data:N*M matrix,N - hotspots number,M - sample num;
label:M*1: The category of the M samples (1,2,3,etc)
p：N*1：p-value of the hotspots in anvoa test
optional:col:X*1,cell of string, describes the color of the X types, if not exist, the color will be randomly assigned.

distance_type:1, clustergram using spearman Correlation coefficient as similarity (distance evalulation method) and 'weighted' as linkage algorithm;
2,clustergram using euclidean as similarity (distance evalulation method) and 'ward' as linkage algorithm

p_cut:the cut-off to choose the features (i.e. 0.01).
feature_num:The number of the type-specific hotspots in each type.
e.g.cluster_analysis_ANVOA(1,0.01,5000,'ANVOA_IFS_matrix.mat'), indicates based on the anvoa filtered matrix, p-value <= 0.01 will be used, and for each category of the data,
5000 hotspots will be selected by the z-score difference, based on these hotspots, clustergram, with spearman Correlation coefficient as similarity (distance evalulation method) and 'weighted' as linkage algorithm, will be used to do hierarchical clustering
After running this comand, in the current workplace, there would two figures (clustering_result.fig and clustering_result.pdf) and one .mat file contain the matrix that was used to do clustering analysis and the category information of these samples (cluster_anvoa_data.mat).

