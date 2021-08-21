function CRAG(data_name,peak_type,varargin)
%%%%%%Call hotspots based on fragment information saved in txt file in folder 'data_name' (i.e. BH01).
%%%%%%Call hotspots based on coverage and fragment length. (Peak_type==1: IFS. Peak_type==2: GC bias corrected IFS.)
%%%optional parameters:
%%fdr: cut-off. default: 0.20
%%distance: Distance cut-off to merge the significant regions nearby. default:200
%%enrichment:whether or not do enrichment for the hotspots:  argument: 'enrichment', 0 or 1. default: 1.

%%%%%%Add the path of the funtions used in this pipeline as workplace
current_path=pwd;
lo=strfind(current_path,'/');
parent_path=current_path(1,1:(lo(end)-1));
addpath(genpath(parent_path));


if ischar(peak_type)
    peak_type=str2double(peak_type);
end
if nargin<2;error('There should be at least two input parameters');end
if nargin==3;error('Input parameters error£¡');end

fdr=0.2;
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


%%Use function 'txt_read' to read fragment informaion from txt file
%%and map the fragment information. All the data was saved in 'data_n/' as 'chri.mat'
parfor i=1:22
    txt_read(data_name,i);
end


%%The output was saved in 'outpeak_n/' as 'chri_peak.mat'.
parfor i=1:22
    Hotspot_call(data_name,peak_type,i,fdr);
end

%%%Merge all the hotspots and the output was saved in 'result_n/' as 'peak_all.mat'
Hotspot_merge(data_name,distance);

%%Write the hotspots in 'result_n/' as 'hotspots.bed'
Hotspot_out(data_name);

%%Plot the fragment pattern around hotspots. The output was saved in
%%'result_n/' as 'peak_plot.fig/pdf'
Hotspot_IFS_plot(data_name);

%%%%%Do enrichment for the hotspots
if (enrichment==1)
    %%Plot the hotspot pattern around TSS. The output was saved in
    %%'result_n/' as 'TSS_plot.fig/pdf'
    Hotspot_TSS(data_name);
    
    %%Plot the hotspot pattern around CTCF. The output was saved in
    %%'result_n/' as 'CTCF_plot.fig/pdf'
    Hotspot_CTCF(data_name);
    
    %%%%%%%%%Calculate the overlap number of hotspots for each chromHMM, which
    %%%%%%%would laters used to chromHMM enrichment
    %%%%%The output was saved in 'result_n/' as 'enri.mat'
    Hotspot_enrichment(data_name);
    
    %%%%%%%%%Calculate the overlap number of random hotspots for each chromHMM, which
    %%%%%%%would laters used to chromHMM enrichment
    %%%%%The output was saved in 'result_n/' as 'enri_randm.mat'
    Hotspot_rand_enrichment(data_name);
    
    %%%%%%%Use fisher exact test to test do the chromHMM enrichment, based on
    %%%%%%%the overlap value of the hotspots with each chromHMM, and the overlap value of the random hotspots with each chromHMM
    %%%%%%%The output was saved in 'result_n/' as 'Enrichment_plot.fig'
    Hotspot_enrichment_fisher(data_name);
end

end
