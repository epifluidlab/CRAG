function IFS_multiClass_classification(input_file_path,fold_number,class_num)
%%%%%%%Do multi-classification classification for IFS matrix using liner svm

%input_file_path, the file name saving the IFS matrix for all the fold
%fold_number: The number of the k-fold validation
%class_num: The types of samples in this classification


if ischar(fold_number)
    fold_number=str2double(fold_number);
end

if ischar(class_num)
    class_num=str2double(class_num);
end

out_file_name=strcat(input_file_path,'/classfication_result.mat');
pre_score=[];

for i=1:fold_number
     path_in=strcat(input_file_path,'/');
    path_in=strcat(path_in,num2str(i));
    %%%%%%%load test data
    test_file=strcat(path_in,'/result_n/test_norm_data.mat');
    load (test_file);
    test_data=ma;
    
    
    b_str='=ma;';
    a_str='train_data';
    
    clear pattern;
    
    %%%%%load the matrix for each type of samples in the training set and
    %%%%%find the sample size of each group
    for j=1:class_num
        temp_file=strcat('/result_n/train_norm_data',num2str(j));
        data_in=strcat(path_in,temp_file); %%The matrix for jth class
        load (data_in);
        t_str=strcat(a_str,num2str(j));
        t_str=strcat(t_str,b_str);
        eval(t_str);
        sample_num(j,1)=length(ma(1,:));
    end
    
    %%%%%%%down-sample all the froup to the sample size of the smallest
    %%%%%%%group to calculate the centroid of each group
    
    clear ma;
    
    sample_size=min(sample_num);
    
     %%%%%%%Calculate the centroid of each group.
    
    for j=1:class_num
        temp_file=strcat('/result_n/train_norm_data',num2str(j));
        data_in=strcat(path_in,temp_file); %%The matrix for jth class
        load (data_in);
        t_str=strcat(a_str,num2str(j));
        t_str=strcat(t_str,b_str);
        eval(t_str);
        
        c_num=length(ma(1,:));
        if (c_num>sample_size)
            loc=randperm(c_num)';
            ma=ma(:,loc(1:sample_size,1));
        end
        pattern(:,j)=mean(ma,2);
        clear ma;
    end
    
   
    %%%%Top 2 predictor
    m=length(test_data(1,:));
    clear lab;
    for j=1:m
        clear temp;
        for w=1:class_num
            co=corr(pattern(:,w),test_data(:,j),'Type','Spearman');
            temp(w,1)=co;
        end
        [~,in]=sortrows(temp,-1);
        lab(j,1)=in(1,1); %%The top 2 predictors for the sample
        lab(j,2)=in(2,1);
    end
    
    %%%%%%%%Binary classifier for top 2 predictor
    train_data=[];
    a=0;
    clear train_label;
    
    for j=1:class_num
        t_str=strcat(a_str,num2str(j));
        temp_comand_A='temp_data=';
        temp_comand_B=strcat(t_str,';');
        temp_comand=strcat(temp_comand_A,temp_comand_B);
        eval(temp_comand);
        train_data=[train_data temp_data];
        b=a+length(temp_data(1,:));
        train_label((a+1):b,1)=j;
        a=b;
    end
    
    
    %%%%%%Train all the pairwise binary classification model (tree)
    for j=1:(class_num-1)
        for k=(j+1):class_num
            X=[train_data(:,train_label==j) train_data(:,train_label==k)]';
            Y=[j*ones(sum(train_label==j),1);k*ones(sum(train_label==k),1)];
            nb(j,k).mb=fitctree(X,Y);
            nb(k,j).mb=nb(j,k).mb;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m=length(test_data(1,:));
    test_data=test_data';
    for j=1:m
        label=predict(nb(lab(j,1),lab(j,2)).mb,test_data(j,:));
        %%%%%%%%%Predict top 1 possible class
        lab(j,3)=label;
    end
    pre_score=[pre_score;lab test_label i*ones(m,1)];
end

%%%%%%Output the prediction result
%%The furst two columns: top 2 predictor predicted by centroid distance
%%The third column: The top 1 predictor predicted by decision tree
%%The forth column: The true class label
%%The fifth column: The fold id
Readme='The furst two columns: top 2 predictor predicted by centroid distance; The third column: The top 1 predictor predicted by decision tree.The forth column: The true class label;The fifth column: The fold id';
save((out_file_name),'pre_score','Readme');

end
