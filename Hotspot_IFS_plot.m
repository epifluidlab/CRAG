function Hotspot_IFS_plot(data_name)
%%%%%%%%%%Plot the fragment pattern around hotspots
%%%%%%%%%%data_name: file name of the data set, i.e. BH01

file_path = strcat('./',data_name);
out_res=strcat(file_path,'/result_n/');
file_path = strcat(file_path,'/data_n/');

listing = dir(file_path);
path(path,file_path);

cou=0;
for i=1:length(listing)   %%obtain all the names of the usefule files
    lo=strfind(listing(i).name,'ch');
    if ~isempty(lo)
        cou=cou+1;
        name{cou,1}=listing(i).name;
    end
end


peak_file=strcat(out_res,'peak_all.mat');
load (peak_file);

num=length(name);

res=zeros(4001,3);
res(:,1)=(1:4001)';
cou=0;
for ii=1:num  %%%%%Get the fragment information of all the chromosomes
    
    te=name{ii,1};
    lo=strfind(te,'.');
    chr=str2double(te(1,4:(lo-1)));
    
    t_name=file_path;
    read_file=strcat(t_name,name{ii,1});
    
    load (read_file); %%Load the fragment information of the current chromosome
    peak_small=peak_a(peak_a(:,1)==chr,2:3);  %%The hotspots regions
    
    
    if isempty(peak_small)
        n=0;
    else
        n=length(peak_small(:,1));
    end
    
    for i=1:n
        cou=cou+1;
        base=ceil((peak_small(i,1)+peak_small(i,2))/2);   %%midpoint of the current hotspot
        left=max(base-2000,1);   %%-2000bp
        right=min(base+2000,region_len);   %%+2000bp
        num1=base-left;
        num2=right-base;
        %%Obtain all the IFS information around hotspots (-2000bp - 2000bp)
        res((2001-num1):(2001+num2),2)=res((2001-num1):(2001+num2),2)+loc(left:right,1);
        res((2001-num1):(2001+num2),3)=res((2001-num1):(2001+num2),3)+1;
    end
end


title_name='Distribution of IFS with different distances to hotspots (Whole Genome)';
y_name='IFS';

res(:,2)=res(:,2)./res(:,3);  %% For each base, calculate the average coverage (length)

control_res=[0 min(res(:,2));0 max(max(res(:,2)))];   %%Line in the midpoint of the hotspots

plot(res(:,1)-2001,res(:,2),'g-',control_res(:,1),control_res(:,2),'r-'); %%Plot the fragment pattern around hotspots
xlabel('Hotspot center');
ylabel(y_name);
title(title_name);

h=gcf;
name_figure=strcat('IFS_plot','.fig');   %%Save the figure
name_figure=strcat(out_res,name_figure);
saveas(h,name_figure);

name_figure2=strcat('IFS_plot','.pdf');   %%Save the figure
name_figure2=strcat(out_res,name_figure2);
saveas(h,name_figure2);


close(gcf);

end

