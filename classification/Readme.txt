(1)The folder contains all the code for classification (two-category classification and TOO).
(2)To run this pipeline, we should run the script 'bam_read.py' to read the fragment for all the samples and write it to txt files.  
For example:python bam_read.py -in C309.filter.bam -out PNAS/C309 could be used to write all the fragment information into the path PNAS/C309. We recommend put all the samples from the same data set into the same folder. For example, when we analysis the Sun et. al. data set, we should put the fragment information of all the samples in the folder PNAS.
(3)After that, IFS_write.m could be used to get the IFS socre for all the chromosomes of the sample. 
For example:IFS_write('PNAS/C309').
(4)k_fold_call.m could be used to call hotspots for each category in a specific fold in k-fold cross validation.
In the function (k_fold_call(fold_info,input_path,out_file_name,peak_type,ith,varargin)), the parameters are shown as follow:

fold_info: an file name of excel (.xlsx or xls) file, which contain all the sample names (string),category (string) and fold id (interger,i.e.1,2,3...), or a cell variable (N*3) contain all this information. For example, 'HCC_healthy_fold.xlsx', which contain three column, in each row,the first column is its sample name, the second column is its category ('HCC' or 'healthy') and the last column is the fold id it begongs to (i.e. 5).
input_path: the path that the fragment information of these samples located in (i.e.PNAS).
out_name_name: The folder the all the output should be saved (i.e. HCC_healthy).
peak_type: 1-call hotspots without GC bias correction; 2- call hotspots based on IFS after GC bias correction
ith: call hotspots for which fold in the cross validation (i.e. 2).

optional parameters (the same with hotspot calling in single sample):
global_p: global p-value cut off
local_p: p-value cut-off for local test
fdr: cut-off
distance: Distance cut-off to merge the significant regions nearby.
enrichment:whether or not do enrichment for the hotspots: 

An example to run the pipeline is: k_fold_call('HCC_healthy_fold.xlsx','PNAS','HCC_healthy',1,2,'enrichment',0);
After run the function, in the folder of 'HCC_healthy', the program will produce a file 'HCC_healthy/2/HCC/, and in HCC_healthy/2/HCC/result_n/, four files would exist (peak_all.mat,peaks.bed) and two files of the fragmentation pattern around hotspots (IFS_plot.fig, IFS_plot.pdf)). Also, the program will produce a file 'HCC_healthy/5/healthy/, and in HCC_healthy/5/healthy/result_n/, four files would exist (peak_all.mat,peaks.bed) and two files of the fragmentation pattern around hotspots (IFS_plot.fig, IFS_plot.pdf))

If do 10-fold validation, the program should be executed 10 times, change ith from 1 to 10.

(4)Function: k_fold_IFS_matrix(fold_info,input_path,out_file_name,peak_type,ith,pos_class) could be used to get the hotspots' IFS of each sample and merge all the samples of a specific category in the training set as a matrix and output the matrix, also the program will output the matrix for the testing samples.
The parameters are shown as follow:

fold_info: an file name of excel (.xlsx or xls) file, which contain all the sample names (string),category (string) and fold id (interger,i.e.1,2,3...), or a cell variable (N*3) contain all this information. For example, 'HCC_healthy_fold.xlsx', which contain three column, in each row,the first column is its sample name, the second column is its category ('HCC' or 'healthy') and the last column is the fold id it begongs to (i.e. 5).
input_path: the path that the fragment information of these samples located in (i.e.PNAS).
out_name_name: The folder the all the output should be saved (i.e. HCC_healthy).
peak_type: 1-get IFS of the hotspots without GC bias correction; 2- get IFS of hotspots based on IFS after GC bias correctionand correct the IFS.
ith: ontain the IFS matrix for which fold in the cross validation (i.e. 2).
pos_class:a string shows the positive class in the classification, i.e. 'HCC', if this is a fold-validation for multi-classification,spefiific it as 'None'.

e.g. After the command k_fold_IFS_matrix('HCC_healthy_fold.xlsx','PNAS','HCC_healthy',1,2,'HCC'), in the folder of 'HCC_healthy', the program will produce a file 'HCC_healthy/2/result_n/, athree files would exist (train_norm_data1.mat,train_norm_data2.mat,test_norm_data.mat).
In train_norm_data1.mat, there is a M*N matrix, N is the number of HCC samples in the second training set. M is the total hotspot number of HCC samples and healthy samples in this fold.
In train_norm_data2.mat, there is a M*N2 matrix, N2 is the number of Healthy samples in the second training set. M is the total hotspot number of HCC samples and healthy samples in this fold.
In test_norm_data.mat, there is a M*N3 matrix, N3 is the number of samples in the second test data set. M is the total hotspot number of HCC samples and healthy samples in this fold (training set). Also, there is a variable 'test_label' (N3 *1), is the class labels of all the test samples.
If do 10-fold validation, the program should be executed 10 times, change ith from 1 to 10.


If for multi-classification, for example, 'Breast cancer vs. lung cancer vs. ovarian cancer vs. colon cancer vs. Pancreatic cancer', we only need to change a parameter in the function.
e.g.
k_fold_call('multi_cancer_fold.xlsx','multi_cancer','TOO',2,2,'enrichment',0);
k_fold_IFS_matrix('multi_cancer_fold.xlsx','multi_cancer','TOO',2,2,'None').
 
In multi_cancer_fold.xlsx: The first column is the sample names for all the five cancer types, the second column is the category and the third column is the fold id.
the folder 'multi_cancer' should contain all the fragment information files of the samples in these data set.
After run this two command, in TOO/2/Breast cancer/, TOO/2/lung cancer/,TOO/2/ovarian cancer/ and TOO/2/Pancreatic cancer/. there should be a file result_n in each of these folers and all these files should contain the hotspot files which are the same with above. Also, in TOO/2/result_n/, there should exist six matrix files:
train_norm_data1.mat, train_norm_data2.mat, train_norm_data3.mat,train_norm_data4.mat,train_norm_data5.mat,test_norm_data.mat.These matrixs are the IFS data for the five cancer types in the training set and all the samples in the test set.

(5)Function: IFS_binaryClass_classification(input_file_path,fold_number,feature_num) could be used to do classification for two-class problem (liner svm model).

the parameters are shown as follow:
input_file_path: The file of the matrix file, e.g. 'HCC_healthy'
fold_number: How many folds are in the validation, e.g. 10
feature_num: The number of features used.

e.g.
IFS_binaryClass_classification('HCC_healthy',10,30000), after this command, the program would perform 10-fold validation classification in the folder 'HCC_healthy', using the most stable 30,000 hostpots as features, and liner svm as model. At last, in the floder 'HCC_healthy/', there would be two figures (ROC.fig and ROC.PDF), also there would be a .mat file (classfication_result.mat), which contains all the detail information in the classification.

(6) Function: IFS_multiClass_classification(input_file_path,fold_number,class_num,feature_num) could be used to do multi-classification for multi-classification problem (TOO).

the parameters are shown as follow:
input_file_path: The file of the matrix file, e.g. 'TOO'
fold_number: How many folds are in the validation, e.g. 10
class_num: Howm many categories are in this classification, e.g. 5
feature_num: The number of features used (decision tree).

e.g.
IFS_multiClass_classification('TOO',10,5,100000), after this command, the program would perform 10-fold validation classification in the folder 'TOO', using centroid distance to do top 2 predictions, and then the most stable 100,000 hostpots as features, and decision tree as model were used to decide top 1 prediction. At last, in the floder 'TOO/', there would be a .mat file (classfication_result.mat), which contains all the detail information in the classification:
%%The furst two columns: top 2 predictor predicted by centroid distance
%%The third column: The top 1 predictor predicted by decision tree
%%The forth column: The true class label
%%The fifth column: The fold id


