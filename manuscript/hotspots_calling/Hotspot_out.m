function Hotspot_out(data_name)

file_path = strcat('./',data_name);

out_res= strcat(file_path,'/result_n/');

te1=strcat(out_res,'/peak_all.mat');
load (te1);  %%load hotspot file

out_path=strcat(out_res,'peaks.bed');

fp = fopen(out_path,'wt');

peak_a=sortrows(peak_a,1);

n=length(peak_a);

fprintf(fp, '%s	', 'Chromosome id');
fprintf(fp, '%s	', 'Start'); %data in bed file was zero based
fprintf(fp, '%s	', 'End');
fprintf(fp, '%s	', 'IFS');
fprintf(fp, '%s	', 'Global p-value');
fprintf(fp, '%s	', 'Local p-value');
fprintf(fp, '%s\n', 'FDR');

for i=1:n
    chrm=strcat('chr',num2str(peak_a(i,1)));
    fprintf(fp, '%s	', chrm);
    fprintf(fp, '%d	', peak_a(i,2)-1); %data in bed file was zero based
    fprintf(fp, '%d	', peak_a(i,3));
    fprintf(fp, '%d	', peak_a(i,4));
    fprintf(fp, '%d	', peak_a(i,5));
    fprintf(fp, '%d	', peak_a(i,6));
    fprintf(fp, '%d\n', peak_a(i,7));
end
fclose(fp);
end





