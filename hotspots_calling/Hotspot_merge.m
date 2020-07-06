function Hotspot_merge(data_name,window_size)
%%Merge the peaks nearby from all the chromosomes
if (nargin==1)
    window_size=200;   %%threshold to merge the peaks.
end

file_path = strcat('./',data_name);
out_res=strcat(file_path,'/result_n/');
file_path = strcat(file_path,'/outpeak_n/');
listing = dir(file_path);
path(path,file_path);


system(['mkdir ' out_res]);


cou=0;
for i=1:length(listing)   %%obtain all the names of the usefule files
    lo=strfind(listing(i).name,'ch');
    if ~isempty(lo)
        cou=cou+1;
        name{cou,1}=listing(i).name;
    end
end


n=length(name);
cou1=0;
for i=1:n
    temp=name{i,1};
    lo2=strfind(name{i,1},'_');
    lo1=strfind(name{i,1},'chr');
    id=str2double(temp(1,(lo1+3):(lo2-1)));
    
    te1=strcat(file_path,name{i,1});
    load (te1);
    
    if ~isempty(peak)
        final_region=peak;
        num=length(final_region(:,1));
        peak=zeros(num,6);
        k=1;
        count=0;
        while (k<=num)
            w=k;
            
            %%If the distance of two peaks is less than the threshold (i.e.200bp)
            while ((k+1) <= num) && ((final_region(k,2)+window_size) > final_region(k+1,1))
                k=k+1;
            end
            count=count+1;
            peak(count,1)= final_region(w,1);
            peak(count,2)= final_region(k,2);
            peak(count,3)=mean(final_region(w:k,3));  %%average the IFS in all the windows
            peak(count,4)= max(final_region(w:k,4));  %%Taking the max global p-value in all the region as teh final p-value
            peak(count,5)= max(final_region(w:k,5));  %%Taking the max local p-value in all the region as teh final p-value
            peak(count,6)= max(final_region(w:k,6));  %%Taking the max q-value in all the region as teh final q-value
            k=k+1;
        end
        
        if (count == 0) || (peak(count,2) ~= final_region(num,2))
            count=count+1;
            peak(count,1:6)= final_region(num,1:6);
        end
        
        if (count < num)
            peak((count+1):num,:)=[];
        end
    end
    
    
    num1=length(peak(:,1));
    peak_a((cou1+1):(cou1+num1),1)=id; %%chromosome id
    peak_a((cou1+1):(cou1+num1),2:7)=peak(:,1:6);
    cou1=cou1+num1;
end

peak_a=sortrows(peak_a,1);

%%Save the final hotspots
out_res=strcat(out_res,'/peak_all.mat');
save(out_res,'peak_a','-v7.3');

end


