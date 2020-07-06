function [ratio,over_peak]=overlap(peak_file_a,peak_file_b)
%%investigate the ratio of hotspots in peak_file_a that were overlapped
%%by hotspots in peak_file_b.
%%ratio: hotspots in peak_file_a was overlapped.
%%over_peak: the overlapped hotspots

load (peak_file_a);%%load the hotspot file A
BH01=peak_a;
index=(1:length(BH01(:,1)))';
flag=zeros(length(BH01(:,1)),1);

load (peak_file_b); %%load the hotspot file B
PNAS=peak_a;

cou=0;
for i=1:22
    a=BH01(BH01(:,1)==i,2:3);
    t_index=index(BH01(:,1)==i,1);
    
    b=PNAS(PNAS(:,1)==i,2:3);
    [a,in]=sortrows(a,1);
    t_index=t_index(in,1);
    
    
    b=sortrows(b,1);
    
    n=length(a);
    m=length(b);
    
    for j=1:n
        t1=a(j,1);
        t2=a(j,2);
        for k=1:m
            if (t1 >= b(k,1) && t1 <= b(k,2)) || (t2 >= b(k,1) && t2 <= b(k,2)) || (b(k,1) >= t1 && b(k,1) <= t2) || (b(k,2) >= t1 && b(k,2) <= t2)
                %%%%%%%Overlapped
                cou=cou+1;
                flag(t_index(j,1),1)=1;
                break;
            end
            if t2 < b(k,1)
                break;
            end
        end
    end
end

ratio=cou/length(BH01); %%The ratio of hotspots of file A was overlapped
index=index(flag==1,1);
over_peak=BH01(index,:); %%The overlapped hotspots

end

