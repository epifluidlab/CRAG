function Hotspot_call_multi_sample(file_list,input_path,out_name,peak_type,varargin)
%%%%%%%%Call hotspots using all the files in file_list...
%%%%%and save the hotspots in out_name/result_n
%file_list:the name of an excel file (include the Suffix name). Or a
%variable include all the the files
%%input_path: The path of all the files
%%out_name: The output file anme
%%peak_type: 1-call hotspots using IFS without GC bias correction
%%peak_type: 2-call hotspots using IFS after GC bias correction

%%%optional parameters:
%%fdr: cut-off
%%distance: Distance cut-off to merge the significant regions nearby.
%%enrichment:whether or not do enrichment for the hotspots:  argument: 'enrichment', 0 or 1. default: 1.

%%%%%%Add the path of the funtions used in this pipeline as workplace
current_path=pwd;
lo=strfind(current_path,'/');
parent_path=current_path(1,1:(lo(end)-1));
addpath(genpath(parent_path));

if ischar(peak_type)
    peak_type=str2double(peak_type);
end
if ischar (file_list) && (contains(file_list,'.xlsx') || contains(file_list,'.xls'))
    %%%%%The input is an excel file
    [~,file_list]=xlsread(file_list);
end

if nargin<4;error('There should be at least four input parameters');end
if nargin==5;error('Input parameters error!');end

fdr=0.20;
distance=200;
enrichment=1;
para_val=varargin;
while length(para_val)>=2
    prop =para_val{1};
    val=para_val{2};
    if ischar(val)
        val=str2double(val);
    end
    para_val=para_val(3:end);
    switch prop
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

if ~exist(out_name,'dir')
    system(['mkdir ' out_name]);
end

%%%%%%%%merge all IFS in the file list and save them in out_name/data_n
parfor i=1:22
    data_merge(input_path,out_name,file_list,i);
end

%%%%%%%Call hotspots
parfor i=1:22
    Hotspot_call(out_name,peak_type,i,global_p,local_p,fdr);
end
%%%%%%merge hotspots
Hotspot_merge(out_name,distance);
Hotspot_out(out_name);
Hotspot_IFS_plot(out_name);

if (enrichment==1)
    Hotspot_TSS(out_name);
    Hotspot_CTCF(out_name);
    Hotspot_enrichment(out_name);
    Hotspot_rand_enrichment(out_name);
    Hotspot_enrichment_fisher(out_name);
end

end
