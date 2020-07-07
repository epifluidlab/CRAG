function IFS_matrix_obtain(file_list,input_path,out_name,peak_type,peak_file,varargin)
%%%%%Obtain the z-score of IFS for all the samples of a specific type and save the IFS matrix in out_name
%file_list:the name of an excel file (include the Suffix name), or a
%variable include all the the files
%%input_path: The path of all the files
%%out_name: The output file name (matrix file name)
%%peak_type£º Peak_type==1: IFS. Peak_type==2: GC bias corrected IFS.
%%peak_file: The sample name for the hotspots to obtain the IFS score.
%%%optional parameters:
%%Several sample names of several hotspot sets
if nargin<5;error('There should be at least two input parameters');end
para_val=varargin;
peak_list{1,1}=peak_file;
if ~isempty(varargin)
    cou=1;
    for i=1:length(para_val)
        cou=cou+1;
        peak_list{cou,1} =para_val{i};
    end
end

if (ischar(peak_type))
    type=str2double(peak_type);
end

n=length(peak_list);
peak=[];
for i=1:n
    peak_loc=strcat(peak_list{i,1},'/result_n/peak_all.mat');
    load (peak_loc);
    peak_num{i,1}=peak_list{i,1};
    peak_num{i,2}=length(peak_a);
    peak=[peak;peak_a];
end


sample=file_list; %%All the samples in the current type
n=length(sample);
feature_num=length(peak);
ma=zeros(feature_num,n);
for i=1:n
    path_name=strcat(input_path,'/');
    path_name=strcat(path_name,sample{i,1});
    [feature_data,peak_origin]=IFS_data_obtain(path_name,peak,peak_type);
    ma(:,i)=feature_data;
end

save((out_name),'ma','peak_num','sample','peak','peak_origin','-v7.3');
%%%%Save the IFS matrix for the current type
end


