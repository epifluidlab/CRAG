function Hotspot_TSS(data_name)
%%%%%%%Plot the hotspot pattern around TSS.
%%%%%%%%%data_name: The file name of the data set. i.e. BH01

file_path = data_name;
out_res= strcat(file_path,'/result_n/');
file_path = strcat(file_path,'/result_n/peak_all.mat');

load (file_path);

load Basic_info/hg19_TSS.mat;   %%The TSS information

load Basic_info/chromosome_info.mat;  %%The chromosome information

count=zeros(4001,2);
count(:,1)=(1:4001)';

seed=unique(peak_a(:,1));
num=length(seed);
total_count=0;

for ii=1:num
    
    chrm=strcat('chr',num2str(seed(ii,1)));
    loc_t=strcmpi(tss(:,1),chrm);
    if sum(loc_t) > 0
        tss_t=tss(loc_t,:);
        tss_t=sortrows(tss_t,2);
        n=length(tss_t);
        tss_loc=zeros(n,1);
        for i=1:n
            if tss_t{i,4}=='+'
                tss_loc(i,1)=str2double(tss_t{i,3});
                tss_kind{i,1}='+';
            else
                tss_loc(i,1)=str2double(tss_t{i,3});
                tss_kind{i,1}='-';
            end
        end
        [tss_loc,index]=sort(tss_loc);
        flag=ones(n,1);
        
        %%%%%%%%%Filter TSS nearby%%%%%%%%%%%%%%%
        for i=1:(n-1)
            for j=(i+1):n
                if (flag(i,1)==0)
                    break;
                end
                if (tss_loc(j,1)-tss_loc(i,1))<=3000
                    flag(j,1)=0;
                else
                    break;
                end
            end
        end
        tss_loc=tss_loc(flag==1,:);
        tss_kind=tss_kind(index,1);
        tss_kind=tss_kind(flag==1,:);
        
        n=length(tss_loc);
        total_count=total_count+n;
        
        peak=peak_a(peak_a(:,1)==seed(ii,1),2:3);
        if isempty(peak)
            m=0;
        else
            m=length(peak(:,1));
        end
        chr_loc=strcmpi(chrm_id,chrm);
        region_len=chrm_len(chr_loc,1);
        
        for i=1:n
            base=tss_loc(i,1); %%the location of the TSS
            left=max(base-2000,0);
            right=min(base+2000,region_len);
            if tss_kind{i,1}=='+'
                for j=1:m %%For the current TSS, search each hotspot to calculate the hotspots which has overlap within [-2000bp 2000bp] of the TSS
                    if ((peak(j,1) >= left) && (peak(j,1) <= right )) || ((peak(j,2) >= left) && (peak(j,2) <= right )) || ((peak(j,1) <= left) && (peak(j,2) >= right))
                        left_new=max(left,peak(j,1));
                        right_new=min(right,peak(j,2));
                        num1=base-left_new;
                        num2=right_new-base;
                        count((2001-num1):(2001+num2),2)=count((2001-num1):(2001+num2),2)+1;
                    end
                end
            else %%%%%%for TSS -, reverse the upstream and downstream of the hotspots
                for j=1:m
                    if ((peak(j,1) >= left) && (peak(j,1) <= right )) || ((peak(j,2) >= left) && (peak(j,2) <= right )) || ((peak(j,1) <= left) && (peak(j,2) >= right))
                        left_new=max(left,peak(j,1));
                        right_new=min(right,peak(j,2));
                        num1=base-left_new;
                        num2=right_new-base;   %%%%Reverse%%%%%%%%%%%%%%%
                        count((4002-2001-num2):(4002-2001+num1),2)=count((4002-2001-num2):(4002-2001+num1),2)+1;
                    end
                end
            end
        end
    end
end

count(:,2)=count(:,2)/total_count;  %%For each base, calculate the fraction of TSS that has overlap with the hotspots

control_res=[0 min(count(:,2));0 max(count(:,2))];
plot(count(:,1)-2001,count(:,2),'g-',control_res(:,1),control_res(:,2),'r-');  %% Plot the pattern around TSS
xlabel('TSS');
ylabel('Fraction of TSS');
title('Distribution of hotspot density with different distances to TSS (Whole Genome)');

h=gcf;
name_figure=strcat('TSS_plot','.fig'); %%save the figures
name_figure=strcat(out_res,name_figure);
saveas(h,name_figure);

name_figure2=strcat('TSS_plot','.pdf'); %%save the figures
name_figure2=strcat(out_res,name_figure2);
saveas(h,name_figure2);

close(gcf);

end
