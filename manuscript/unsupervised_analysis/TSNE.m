function TSNE(distance_type,neighbor_num,p_cut,matrix_file)
%%%%%%%%TSNE visualize the samples from the matrix file.
%%The file should including at least three variable:
%%data:N*M matrix,N - hotspots number,M - sample num;
%%label:M*1: The category of the M samples (1,2,3,etc)
%%p£ºN*1£ºp-value of the hotspots in anvoa test
%%optional:col:X*1,cell of string, describes the color of the X types, if
%%not exist, the color will be randomly assigned.

%%%%%%%%distance_type =1, the distance among the samples will be evalulate by
%%spearman Correlation c oefficient
%%distance_type =2: the distance among the samples will be evalulate by euclidean distance
%%neighbor_num: the number of neighbors

%%p_cut, the cut-off to choose the features (i.e. 0.01)


if ischar(distance_type)
    distance_type=str2double(distance_type);
end
if ischar(neighbor_num)
    neighbor_num=str2double(neighbor_num);
end
if ischar(p_cut)
    p_cut=str2double(p_cut);
end

load (matrix_file);
if iscell(label)
    class_label=str2double(label);
else
    class_label=label;
end
data=data(p<=p_cut,:);

class_num=length(unique(class_label));
if ~exist('col','var')
    coArray={'y','m','c','r','g','b','w','k'}';
    col_loc=randi(8,class_num,1);
    col_loc=[randperm(8)';col_loc];
    col(1:class_num)=coArray(col_loc(1:class_num),1);
end
disp(col);



rng default
%%%%%%%%Using t-sne to get the data of the first two dimensions
if (distance_type==1)
    score = tsne(data','Algorithm','exact','Distance','spearman','Perplexity',neighbor_num);
else
    score = tsne(data','Algorithm','exact','Distance','euclidean','Perplexity',neighbor_num);
end


%%%%%%%%%%%plot the data
for i=1:class_num
    temp=score(class_label==i,:);
    scatter(temp(:,1),temp(:,2),15,col{i},'filled');
    hold on
end
xlabel('tSNE 1');
ylabel('tSNE 2');
h=gcf;
name_figure1=strcat('TSNE_plot','.fig');   %%Save the figure
saveas(h,name_figure1);
name_figure2=strcat('TSNE_plot','.pdf');   %%Save the figure
saveas(h,name_figure2);
%%%%%%Also save the data
save('tsne_score.mat','class_label','score');

end

