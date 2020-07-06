function cluster_analysis_ANVOA(distance_type,p_cut,feature_num,matrix_file)
%%%%%%%%clustering the samples from the matrix file.
%%The file should including at least three variable:
%%data:N*M matrix,N - hotspots number,M - sample num;
%%label:M*1: The category of the M samples (1,2,3,etc)
%%p£ºN*1£ºp-value of the hotspots in anvoa test
%%optional:col:X*1,cell of string, describes the color of the X types, if
%%not exist, the color will be randomly assigned.

%%%%%%%%distance_type =1, clustergram using spearman Correlation coefficient
%%as similarity (distance evalulation method) and 'weighted' as linkage algorithm
%%distance_type =2:clustergram using euclidean as similarity (distance evalulation method) and 'ward' as linkage algorithm
%%p_cut, the cut-off to choose the features (i.e. 0.01)
%%feature_num: The number of the type-specific hotspots in each type.


if ischar(distance_type)
    distance_type=str2double(distance_type);
end
if ischar(feature_num)
    feature_num=str2double(feature_num);
end
if ischar(p_cut)
    p_cut=str2double(p_cut);
end

if nargin<4;error('There should be at least four input parameters');end

load (matrix_file); %%load the matrix
if iscell(label)
    class_label=str2double(label);
else
    class_label=label;
end

ma=data(p<=p_cut,:);

index=(1:length(ma(:,1)))';
loc=[];
class_num=length(unique(class_label));

%%%%%%For the significant IFS, using z-score difference to select the
%%%%%%disease-specific hotspots. top feature_num for each disease
for i=1:class_num
    val=abs(mean(ma(:,class_label==i),2)-mean(ma(:,class_label~=i),2));
    t_index=[index val];
    t_index=sortrows(t_index,-2);
    loc=[loc;t_index(1:feature_num,1)];
end
loc=unique(loc);

%%%%%%If the color doesn't exist, random assign the color and disp the
%%%%%%color in the screen
if ~exist('col','var')
    coArray={'y','m','c','r','g','b','w','k'}';
    col_loc=randi(8,class_num,1);
    col_loc=[randperm(8)';col_loc];
    col(1:class_num)=coArray(col_loc(1:class_num),1);
end
disp(col);


%%%%%%%%%The final data for clustering analysis
data=-ma(loc,:);
clear label;
for i=1:length(class_label)
   label{i,1}=num2str(class_label(i,1)); 
end

color_define = cell(size(label));
for i=1:class_num
    color_define(strcmpi(label,num2str(i)))=col(i);
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
name_figure1=strcat('clustering_anvoa_result','.fig'); %%save the figures
saveas(h,name_figure1);
name_figure2=strcat('clustering_anvoa_result','.pdf'); %%save the figures
saveas(h,name_figure2);
save('cluster_anvoa_data.mat','data','label','-v7.3');
%%%%%%%Output the clustering result as .fig file and .pdf file and the feature matrix, as well as label information .
end
