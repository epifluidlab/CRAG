function ANVOA_calculate(p_cut,out_name,matrix_file1,color1,matrix_file2,color2,varargin)
%%%%%%%Merge the IFS score for all the samples and use one-way anvoa
%%%%%%%to select the features
%%p_cut: The p-value cut off in one-way anvoa to select the features
%%out_name: The name of the output matrix (selcted features),class_label of
%%these samples,Color for each category (for future analysis)
%%matrix_file1£ºThe .mat file for the IFS matrix of the first type
%%matrix_file2£ºThe .mat file for the IFS matrix of the second type
%%Options:more the matrix files and the coressponding color

if ischar(p_cut)
    p_cut=str2double(p_cut);
end

if nargin<6;error('There should be at least four input parameters');end

file_list{1,1}=matrix_file1;
file_list{2,1}=matrix_file2;
col{1,1}=color1;
col{2,1}=color2;
cou=2;
para_val=varargin;
while length(para_val)>=2
    cou=cou+1;
    file_list{cou,1} =para_val{1};
    col{cou,1}=para_val{2};
    para_val=para_val(3:end);
end


class_num=cou;
data=[];
sample_num=0;
for i=1:class_num
    load (file_list{i,1});
    data=[data ma];
    for j=1:length(ma(1,:))
        sample_num=sample_num+1;
        label{sample_num,1}=num2str(i);
    end
end


n=length(data(:,1));
p=zeros(n,1);
for i=1:n
  p(i,1) = anova1(data(i,:),label,'off');  
end
data=data(p<=p_cut,:);
p=p(p<=p_cut,1);
%%%%%Using anvoa to select the significant IFS and save the significant
%%%%%data
save((out_name),'data','label','p','col','-v7.3');

end
