function [ratio,over_peak]=overlap_ratio_repeat(peak_file_A,repeat_file_B)
%%%%%This file is for figure 4e
%%investigate the ratio of hotspots in peak_file_a that were overlapped
%%by the regions in repeat file.
%%ratio: hotspots in peak_file_a was overlapped.
%%over_peak: the overlapped hotspots


load (peak_file_A); %%load hotspot file
BH01=peak_a; 
   
load (repeat_file_B); %% load repeat file

PNAS=loca;  %%loca is the regions of the repeat.
over_peak=[];
non_over_peak=[];

for i=1:22
    a=BH01(BH01(:,1)==i,:);
    
    chrm=strcat('chr',num2str(i));
    lo_b=strcmpi(chr_id,chrm);
    b=PNAS(lo_b,:);
    
    b(:,1)=b(:,1)+1; %%%repeat is 0-based
    
    a=sortrows(a,2);
    
    flag=zeros(length(a),1);
    
    b=[ones(length(b),1)*i,b];
    b=sortrows(b,2);
    
    n=length(a);
    m=length(b);
    
    for j=1:n
        t1=a(j,2);
        t2=a(j,3);
        for k=1:m
            if (t1 >= b(k,2) && t1 <= b(k,3)) || (t2 >= b(k,2) && t2 <= b(k,3)) || (b(k,2) >= t1 && b(k,2) <= t2) || (b(k,3) >= t1 && b(k,3) <= t2)
                %%overlapped hotspots
                flag(j,1)=1;
                break;
            end
            if t2 < b(k,2)
                break;
            end
        end
    end
    over_peak=[over_peak;a(flag==1,:)];
    non_over_peak=[non_over_peak;a(flag==0,:)]; 
end

ratio=length(over_peak)/length(BH01);

end
