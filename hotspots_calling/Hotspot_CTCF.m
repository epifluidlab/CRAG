function Hotspot_CTCF(data_name)
%%%%%%%Plot the hotspot pattern around CTCF.
%%%%%%%%%data_name: The file name of the data set. i.e. BH01


file_path = data_name;
out_res= strcat(file_path,'/result_n/');
file_path = strcat(file_path,'/result_n/peak_all.mat');

load (file_path);

load Basic_info/hg19_CTCF.mat;  %%The CTCF information

load Basic_info/chromosome_info.mat; %%The chromosome information

count_s=zeros(4001,2);
count_s(:,1)=(1:4001)';

seed=unique(peak_a(:,1));
num=length(seed);

total_count=0;

for ii=1:num
    chrm=strcat('chr',num2str(seed(ii,1)));
    
    loc_t=strcmpi(CTCF(:,1),chrm);
    chr_loc=strcmpi(chrm_id(:,1),chrm);
    region_len=chrm_len(chr_loc,1);
    
    
    if sum(loc_t) > 0
        CTCF_t=CTCF(loc_t,:);
        CTCF_t=sortrows(CTCF_t,2);
        n=length(CTCF_t);
        CTCF_loc=zeros(n,1);
        for i=1:n
            if CTCF_t{i,4}=='+'
                CTCF_loc(i,1)=fix((str2double(CTCF_t{i,2})+1+str2double(CTCF_t{i,3}))/2);  %% bed file is zeo-based. Calculate the midpoint of the CTCF.
                CTCF_kind{i,1}='+';
            else
                CTCF_loc(i,1)=fix((str2double(CTCF_t{i,2})+1+str2double(CTCF_t{i,3}))/2);  %% bed file is zeo-based. Calculate the midpoint of the CTCF.
                CTCF_kind{i,1}='-';
            end
        end
        [CTCF_loc,index]=sort(CTCF_loc);
        flag=ones(n,1);
        
        %%%%%%%%%%%%%%%%Filter CTCF nearby%%%%%%%%%%%%%%%%%
        for i=1:(n-1)
            for j=(i+1):n
                if (flag(i,1)==0)
                    break;
                end
                if (CTCF_loc(j,1)-CTCF_loc(i,1))<=3000
                    flag(j,1)=0;
                else
                    break;
                end
            end
        end
        CTCF_loc=CTCF_loc(flag==1,:);
        
        CTCF_kind=CTCF_kind(index,1);
        CTCF_kind=CTCF_kind(flag==1,:);
        
        
        n=length(CTCF_loc);
        total_count=total_count+n;
        
        peak=peak_a(peak_a(:,1)==seed(ii,1),2:3);
        peak_small=peak;
        if isempty(peak_small)
            m=0;
        else
            m=length(peak_small(:,1));
        end
        for i=1:n
            base=CTCF_loc(i,1);
            left=max(base-2000,0);
            right=min(base+2000,region_len);
            if CTCF_kind{i,1}=='+'
                for j=1:m %%For the current CTCF, search each hotspot to calculate the hotspots which has overlap within [-2000bp 2000bp] of the TSS
                    if ((peak_small(j,1) >= left) && (peak_small(j,1) <= right )) || ((peak_small(j,2) >= left) && (peak_small(j,2) <= right )) || ((peak_small(j,1) <= left) && (peak_small(j,2) >= right ))
                        left_new=max(left,peak_small(j,1));
                        right_new=min(right,peak_small(j,2));
                        num1=base-left_new;
                        num2=right_new-base;
                        count_s((2001-num1):(2001+num2),2)=count_s((2001-num1):(2001+num2),2)+1;
                    end
                end
            else %%%%%%for CTCF -, reverse the upstream and downstream of the hotspots
                for j=1:m
                    if ((peak_small(j,1) >= left) && (peak_small(j,1) <= right )) || ((peak_small(j,2) >= left) && (peak_small(j,2) <= right )) || ((peak_small(j,1) <= left) && (peak_small(j,2) >= right ))
                        left_new=max(left,peak_small(j,1));
                        right_new=min(right,peak_small(j,2));
                        num1=base-left_new;
                        num2=right_new-base;   %%%%%%%Reverse%%%%%%%%%%%%%%%%%%%%%%%%
                        count_s((4002-2001-num2):(4002-2001+num1),2)=count_s((4002-2001-num2):(4002-2001+num1),2)+1;
                    end
                end
            end
        end
    end
end

count=count_s;
count(:,2)=count(:,2)/total_count; %%For each base, calculate the fraction of CTCF that has overlap with the hotspots

res=[0 min(count(:,2));0 max(count(:,2))];
plot(count(:,1)-2001,count(:,2),'g-',res(:,1),res(:,2),'r-');  %% Plot the pattern around CTCF
xlabel('CTCF');
ylabel('Fraction of CTCF');
title('Distribution of Hotspot density with different distances to CTCF motifs (Whole Genome)');

h=gcf;
name_figure=strcat('CTCF_plot','.fig'); %%save the figures
name_figure=strcat(out_res,name_figure);
saveas(h,name_figure);

name_figure2=strcat('CTCF_plot','.pdf'); %%save the figures
name_figure2=strcat(out_res,name_figure2);
saveas(h,name_figure2);

close(gcf);

end
