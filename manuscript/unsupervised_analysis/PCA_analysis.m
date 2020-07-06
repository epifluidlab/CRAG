function PCA_analysis(matrix_file,color,varargin)
%%%%%%%%Do PCA for the samples from different matrix files (each matrix file contain one type of sampples).
%%matrix_file£ºThe .mat file for the IFS matrix of this type
%%color:after clustering, the color do denote this type
%%Options: 2*N (N>=1) variable, the matrix file and its corresponding color

if nargin<2;error('There should be at least two input parameters');end
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
        label(sample_num,1)=i;
    end
end

ma=data';

%%%%%%%%%%%PCA
[coeff,score,latent,tsquared,explained] = pca(ma);

%%%%%%%%%%%plot the data
for i=1:type_num
    temp=score(label==i,:);
    scatter(temp(:,1),temp(:,2),15,col{i},'filled');
    hold on
end

xlabel('PCA1');
ylabel('PCA2');
h=gcf;
name_figure1=strcat('PCA_plot','.fig');   %%Save the figure
saveas(h,name_figure1);
name_figure2=strcat('PCA_plot','.pdf');   %%Save the figure
saveas(h,name_figure2);
%%%%%%Also save the data
save('PCA_score.mat','label','score');



end




