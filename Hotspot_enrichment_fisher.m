function Hotspot_enrichment_fisher(data_name)
%%%%%%%Use fisher exact test to test do the chromHMM enrichment, based on
%%%%%%%the overlap value of the hotspots with each chromHMM, and the overlap value of the random hotspots with each chromHMM
%%%%%%%%%data_name: The file name of the data set. i.e. BH01

file_path = strcat('./',data_name);
out_res= strcat(file_path,'/result_n/');


te1=strcat(out_res,'peak_all.mat');
load (te1); %%load all the hotspots data
peak_num=length(peak_a);  %%Hotspot number

te2=strcat(out_res,'enri.mat');   %% The overlap value of true hotspots
te3=strcat(out_res,'enri_randm.mat');  %%The overlap value of random hotspots
load (te2);
load (te3);

p_ma=ones(15,length(ma(1,:)));

for i=1:length(ma(:,1))
    for j=1:length(ma(1,:))
        if ma(i,j)~=0
            x_table=[ma(i,j) peak_num-ma(i,j);ma_c(i,j) peak_num-ma_c(i,j)];
            [~,p,~] = fishertest(x_table,'tail','right');  %%fisher exact test to test whether the overlap of the hotspot to the current chromHMM is significantly higher than the random hotspots
            if ma_c(i,j)~=0
                p_ma(i,j)=p;
            else
                p_ma(i,j)=1;
            end
        else
            p_ma(i,j)=1;
        end
    end
end

for i=1:length(ma(:,1))
    for j=1:length(ma(1,:))
        if p_ma(i,j)==0 %% The max precision
            p_ma(i,j)=10e-300;
        end
    end
end

[id_index,y_name]=xlsread('./Basic_info/cell.xlsx');   %%Cell types
[~,x_name]=xlsread('./Basic_info/state.xlsx'); %%chromHmm names

p_ma=p_ma(:,id_index);


xvalues = x_name;
yvalues = y_name(:,2);
m=heatmap(yvalues,xvalues,-log10(p_ma),'Colormap',jet);  %%Heatmap for all the enrichment value.

m.Title = 'Enrichment alanysis of our hotspots with different chromatin states';
m.YLabel = 'chromatin states';

h=gcf;
name_figure=strcat('Enrichment_plot','.fig');  %%Save the enrichment result.
name_figure=strcat(out_res,name_figure);
saveas(h,name_figure);

name_figure2=strcat('Enrichment_plot','.pdf');  %%Save the enrichment result.
name_figure2=strcat(out_res,name_figure2);
saveas(h,name_figure2);
close(gcf);



end