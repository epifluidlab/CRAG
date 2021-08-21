function IFS_binaryClass_classification(input_file_path,fold_number)
%%%%%%%Do binary classification for IFS matrix using liner svm
%input_file_path, the file name saving the IFS matrix for all the fold
%fold_number: The number of the k-fold validation

if ischar(fold_number)
    fold_number=str2double(fold_number);
end


out_file_name=strcat(input_file_path,'/classfication_result.mat');

pre_score=[];
auc_score=zeros(fold_number,1);
sensitivity(1,1)=1;
sensitivity(1,2)=0.99;
sensitivity(1,3)=0.98;
sensitivity(1,4)=0.95;
sensitivity(1,5)=0.85;
for i=1:fold_number
    path_in=strcat(input_file_path,'/');
    path_in=strcat(path_in,num2str(i));
    data_in=strcat(path_in,'/result_n/train_norm_data1.mat'); %%The matrix for the positive class
    load (data_in);
    train_data1=ma;
    pos_num=length(ma(1,:));
    clear train_label;
    train_label(1:pos_num,1)=1; %%class label for positive samples
    
    data_in=strcat(path_in,'/result_n/train_norm_data2.mat'); %%The matrix for the negative class
    load (data_in);
    train_data2=ma;
    neg_num=length(ma(1,:));
    train_label((pos_num+1):(pos_num+neg_num),1)=-1; %%class label for negative samples
    
    data_in=strcat(path_in,'/result_n/test_norm_data.mat'); %%The matrix for the test samples
    load (data_in);
    test_data=ma;
    

    
    %%Only keep the IFS for the selected features 
    train_data=[train_data1 train_data2];
    
    
    %%%%%%%%%%%%Divide the x axis to 0,0.01,0.02,...1
    intervals= linspace(0, 1, 101);
    
    
    X=train_data';
    Y=train_label(:,1);
    %%%%%train Liner SVM model
    SVMModel = fitcsvm(X,Y,'ClassNames',[-1,1]);
    test=test_data';
    %%%%%%%%%predict the test data set
    [Label,score] = predict(SVMModel,test);
    pre_score=[pre_score;[Label score(:,2) test_label ones(length(test_label),1)*i]];
    [X,Y,~,AUC] = perfcurve(test_label,score(:,2),'1');
   
    %%%%%%%%Use function to estimate the sensitivity from 1-specificity [0 0.01 0.02...1]
    val=data_point_estimate(X,Y,intervals);
    if i==1
        mean_curve= val/fold_number;
    else
        mean_curve= mean_curve+ val/fold_number;
    end
    sensitivity(i+1,1)=val(1);
    sensitivity(i+1,2)=val(2);
    sensitivity(i+1,3)=val(3);
    sensitivity(i+1,4)=val(6);
    sensitivity(i+1,5)=val(16);
    auc_score(i,1)=AUC;
    
end
X=[0 intervals]';
Y=[0;mean_curve];
%%%%%%%%%%Output the X axis values and Y axis values for ROC
%%%%%%%%%%Output the auc in the K folds
%%%%%%%%%%Output sensitivity in each fold at several speccificities.
%%%%%%%%%%Output The predicted label, predicted score, true label,fold id
plot(X,Y);
xlabel('1 - Specificity');
ylabel('Sensitivity');
h=gcf;
name_figure=strcat('ROC','.fig');   %%Save the figure
out_res=strcat(input_file_path,'/');
name_figure=strcat(out_res,name_figure);
saveas(h,name_figure);

name_figure2=strcat('ROC','.pdf');   %%Save the figure
name_figure2=strcat(out_res,name_figure2);
saveas(h,name_figure2);
Readme='Variable_describe:X and Y are values in x axis and y axis in ROC; auc_score contain all the auc in k folds;in sensitivity, the first row is the specificity,2 - fold_num+1 rows are the sensitivity @specificity of 1st row in k-folds;in pre_score, the first column is the predicted label, the second column is the predicted score, the third column is the true label, the forth column is the fold id.';

save((out_file_name),'X','Y','auc_score','sensitivity','pre_score','Readme');


end
