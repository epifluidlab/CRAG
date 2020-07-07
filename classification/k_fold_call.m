function k_fold_call(fold_info,input_path,out_file_name,peak_type,ith,varargin)
%%%%%%%%Call hotspots for the hotspots in the training set of ith fold
%%%%%and save the hotspots in out_name/ith/diseasetype/result_n
%fold_info:the name of an excel file (include the Suffix name). Or a
%variable include all the the files
%the excel file or the variable should include three columns:
%1st, sample_name (String: e.g. C302); 2nd, disease type (String: e.g. breast); 3rd, an interger indicate which fold the sample is belonged to.

%%input_path: The path of all the sample files

%%out_file_name: The folder the all these files should be saved.
%%peak_type: 1-call hotspots using IFS without GC bias correction
%%peak_type: 2-call hotspots using IFS after GC bias correction

%%ith: the id of the current fold.


%%%optional parameters:

%%global_p: global p-value cut off
%%local_p: p-value cut-off for local test
%%fdr: cut-off
%%distance: Distance cut-off to merge the significant regions nearby.
%%enrichment:whether or not do enrichment for the hotspots:  argument: 'enrichment', 0 or 1. default: 0


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

if nargin<5;error('There should be at least five input parameters');end
if nargin==6;error('Input parameters error!');end

global_p=0.00001;
local_p=0.00001;
fdr=0.01;
distance=200;
enrichment=0;
para_val=varargin;
while length(para_val)>=2
    prop =para_val{1};
    val=para_val{2};
    if ischar(val)
        val=str2double(val);
    end
    para_val=para_val(3:end);
    switch prop
        case 'global_p'
            global_p=val;
        case 'local_p'
            local_p=val;
        case 'fdr'
            fdr=val;
        case 'distance'
            distance=val;
        case 'enrichment'
            enrichment=val;
        otherwise
            error('Input parameters error')
    end
end
path_base=strcat(out_file_name,'/');
out_name=strcat(path_base,num2str(ith));
system(['mkdir ' out_name]);

train_sample=sample_info(fold_id~=ith,:);

type_name=unique(train_sample(:,2));

n=length(type_name(:,1));
out_name=strcat(out_name,'/');
for i=1:n
   dis_out_name{i,1}=strcat(out_name,type_name{i,1});
   system(['mkdir ' dis_out_name{i,1}]);
   loc=strcmpi(train_sample(:,2),type_name{i,1});
   current_sample_list=train_sample(loc,1);
   Hotspot_call_multi_sample(current_sample_list,input_path,dis_out_name{i,1},peak_type,'enrichment',enrichment);
end

end
