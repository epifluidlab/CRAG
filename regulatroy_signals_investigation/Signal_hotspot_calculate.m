function Signal_hotspot_calculate(data_name,signal_type,hotspot_file)
%%%Calculate the histone siganls around hotspots
%%data_name:the sample name of signals, i.e.E29
%%signal_type: the type of histone siganls:DNase; H3K4me1;H3K4me3...
%%H3K9me3;H3K27ac;H3K27me3;H3K36me3
%%hotspot_file:the name of the hotspot sample:i.e.BH01


file_path = strcat('./Basic_info/Histone_Signal/',data_name);
out_res=strcat(file_path,'/');

file_path = strcat(file_path,'/');
file_path = strcat(file_path,signal_type);

listing = dir(file_path);
path(path,file_path);

cou=0;
for i=1:length(listing)   %%obtain all the names of the chromosome id
    lo=strfind(listing(i).name,'ch');
    if ~isempty(lo)
        cou=cou+1;
        name{cou,1}=listing(i).name;
    end
end

peak_file=strcat(hotspot_file,'/result_n/peak_all');

load (peak_file);

t_file=strcat(signal_type,'_');
t_file=strcat(t_file,hotspot_file);


if (cou~=0)
    num=length(name);
    res=zeros(4001,3);
    res(:,1)=(1:4001)'-2001;
    for ii=1:num %%For each chromosome
        te=name{ii,1};
        lo=strfind(te,'.');
        chr=str2double(te(1,4:(lo-1)));
        
        t_name=strcat(file_path,'/');
        read_file=strcat(t_name,name{ii,1});
        
        load (read_file);
        peak_small=peak_a(peak_a(:,1)==chr,2:3);
        n=length(peak_small);
        for i=1:n
            base=ceil((peak_small(i,1)+peak_small(i,2))/2);
            if ((base-2000)>=1) && ((base+2000) <= length(loc))
                left=base-2000;
                right=base+2000;
                %%%%%%%Calculate the siganls around the hotspots [-2000bp 2000bp]
                res(:,2)=res(:,2)+loc(left:right,1);
                res(:,3)=res(:,3)+1;
            end
        end
    end
    
    x=res(:,1); %%x axis (distance to the hotspot)
    y=res(:,2)./res(:,3); %% the average signals around the hotspots
    name_data=strcat(out_res,t_file);
    save((name_data),'x','y');
end

end

