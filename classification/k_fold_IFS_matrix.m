function k_fold_IFS_matrix(fold_info,input_path,out_file_name,peak_type,ith,pos_class)
%%%%%%%%Based the hotspots in the training set of ith fold,obtain the IFS
%%%%%%%%of all the hotspots in all the samples, including in the training
%%%%%%%%set and test set.
%fold_info:the name of an excel file (include the Suffix name). Or a
%variable include all the the files
%the excel file or the variable should include three columns:1th,
%sample_name (String,i,e,C302),2th,disease type (String:i.e.breast),3th,a interger indicate which
%fold the sample is belonged to.

%%input_path: The path of all the sample files

%%out_file_name: The folder the all these files should be saved.
%%peak_type: 1-call hotspots using IFS without GC bias correction
%%peak_type: 2-call hotspots using IFS after GC bias correction

%%ith: the id of the current fold.
%%pos_class:a string shows the positive class in the
%%disease,i.e.breast,if this is a fold-call for
%%multi-classification,spefiific it as 'None'.

%%%%%%Add the path of the funtions used in this pipeline as workplace
current_path=pwd;
lo=strfind(current_path,'/');
parent_path=current_path(1,1:(lo(end)-1));
addpath(genpath(parent_path));



if ischar(peak_type)
    peak_type=str2double(peak_type);
end
if ischar(ith)
    ith=str2double(ith);
end


if ischar (fold_info) && (contains(fold_info,'.xlsx') || contains(fold_info,'.xls'))
    %%%%%The input is an excel file
    [fold_id,sample_info]=xlsread(fold_info);
else
    fold_id=cell2mat(fold_info(:,3));
    if ischar(fold_id)
        fold_id=str2double(fold_id);
    end
    sample_info=fold_info(:,1:2);
end

path_base=strcat(out_file_name,'/');
out_name=strcat(path_base,num2str(ith));
new_out_name=strcat(out_name,'/result_n/');
system(['mkdir ' new_out_name]);

train_sample=sample_info(fold_id~=ith,:);
test_sample=sample_info(fold_id==ith,:);

type_name=unique(train_sample(:,2));
n=length(type_name(:,1));
out_name=strcat(out_name,'/');
peak=[];
for i=1:n
    peak_list{i,1}=strcat(out_name,type_name{i,1});
    peak_loc=strcat(peak_list{i,1},'/result_n/peak_all.mat');
    load (peak_loc);
    peak=[peak;peak_a];
end
peak_a=peak;
peak_out_name=strcat(new_out_name,'/peak_all.mat');
save((peak_out_name),'peak_a');


if strcmpi(pos_class,'None')
    for i=1:n
        %%%%multi-classification, get and save the matrix for each disease type
        temp_out_name=strcat('train_norm_data',num2str(i));
        result_out=strcat(new_out_name,temp_out_name);
        file_loc=strcmpi(train_sample(:,2),type_name{i,1});
        file_list=train_sample(file_loc,1);
        IFS_matrix_obtain(file_list(:,1),input_path,result_out,peak_type,out_name);
    end

    %%%%%%%%%%get the label for test data set and get the matrix
    test_label=zeros(length(test_sample(:,1)),1);
    for j=1:n
        file_loc=strcmpi(test_sample(:,2),type_name{j,1});
        test_label(file_loc,1)=j;
    end
    temp_out_name='test_norm_data.mat';
    result_out=strcat(new_out_name,temp_out_name);
    IFS_matrix_obtain(test_sample(:,1),input_path,result_out,peak_type,out_name);
    load (result_out);
    %%%add the test_label information to the data
    save((result_out),'ma','peak_num','sample','peak','peak_origin','test_label')
else
    %%%%%%%%binary classification, get and save the matrix for pos-class
    %%%%%%%%and negative class
    %%positive
    temp_out_name='train_norm_data1.mat';
    result_out=strcat(new_out_name,temp_out_name);
    file_loc=strcmpi(train_sample(:,2),pos_class);
    file_list=train_sample(file_loc,1);
    IFS_matrix_obtain(file_list,input_path,result_out,peak_type,out_name);
    %%negative class
    temp_out_name='train_norm_data2.mat';
    result_out=strcat(new_out_name,temp_out_name);
    file_list=train_sample(~file_loc,1);
    IFS_matrix_obtain(file_list,input_path,result_out,peak_type,out_name);

    test_label=zeros(length(test_sample(:,1)),1);
    file_loc=strcmpi(test_sample(:,2),pos_class);
    test_label(file_loc,1)=1;
    test_label(~file_loc,1)=-1;
    temp_out_name='test_norm_data.mat';
    result_out=strcat(new_out_name,temp_out_name);
    IFS_matrix_obtain(test_sample(:,1),input_path,result_out,peak_type,out_name);
    load (result_out);
    %%%add the test_label information to the data
    save((result_out),'ma','peak_num','sample','peak','peak_origin','test_label','-v7.3')
end

end
