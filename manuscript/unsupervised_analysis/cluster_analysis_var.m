function cluster_analysis_var(distance_type,feature_num,matrix_file,color,varargin)
%%%%%%%%clustering the samples from different matrix files (each matrix file contain one type of sampples).
%%%%%%%%distance_type =1, clustergram using spearman Correlation coefficient
%%as similarity (distance evalulation method) and 'weighted' as linkage algorithm
%%distance_type =2:clustergram using euclidean as similarity (distance evalulation method) and 'ward' as linkage algorithm
%%feature_num: The number of the hotspots were used for clustering.
%%matrix_file£ºThe .mat file for the IFS matrix of this type
%%color:after clustering, the color do denote this type
%%Options: 2*N (N>=1) variable, the matrix file and its corresponding color

if ischar(distance_type)
    distance_type=str2double(distance_type);
end
if ischar(feature_num)
    feature_num=str2double(feature_num);
end

if nargin<4;error('There should be at least four input parameters');end
if nargin==5;error('Input parameters error£¡');end

file_list{1,1}=matrix_file;
col{1,1}=color;
cou=1;
para_val=varargin;
while length(para_val)>=2
    cou=cou+1;
    file_list{cou,1} =para_val{1};
    col{cou,1}=para_val{2};
    para_val=para_val(3:end);
end
type_num=cou;
data=[];
sample_num=0;
for i=1:type_num
    load (file_list{i,1});
    data=[data ma];
    for j=1:length(ma(1,:))
        sample_num=sample_num+1;
        label{sample_num,1}=num2str(i);
    end
end


index=(1:length(data(:,1)))';
vari=var(data')';
index(:,2)=vari;

index=sortrows(index,-2); % The hotspots were ranked by the variance across all samples

data=-data(index(1:feature_num,1),:);

color_define = cell(size(label));
for i=1:type_num
    color_define(strcmpi(label,num2str(i)))  = col(i,1);
end

s = struct('Labels',label,'Colors',color_define);
%%%%%%%%%clustergram to do clustering analysis for all the samples
if (distance_type==1)
    cgo_all=clustergram(data,'Linkage','weighted','ColumnPDist','Spearman','RowPDist','Spearman','Colormap',parula);%v2
else
    cgo_all=clustergram(data,'Linkage','ward','Colormap',parula);%v2
end

set(cgo_all,'ColumnLabels',label,'ColumnLabelsColor',s,'LabelsWithMarkers',true);
h=plot(cgo_all);
name_figure1=strcat('clustering_result','.fig'); %%save the figures
saveas(h,name_figure1);
name_figure2=strcat('clustering_result','.pdf'); %%save the figures
saveas(h,name_figure2);

save('cluster_var_data.mat','data','label','-v7.3');
%%%%%%%Output the clustering result as .fig file and .pdf file and the feature matrix, as well as label information .

end
