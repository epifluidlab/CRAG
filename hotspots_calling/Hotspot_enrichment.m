function Hotspot_enrichment(data_name)
%%%%%%%Calculate the overlap number of hotspots for each chromHMM, which
%%%%%%%would laters used to chromHMM enrichment
%%%%%%%%%data_name: The file name of the data set. i.e. BH01

file_path = data_name;
out_res= strcat(file_path,'/result_n/');
file_path = strcat(file_path,'/result_n/peak_all.mat');

load (file_path);   %%load all the hotspots data

current_path=pwd;
lo=strfind(current_path,'/');
parent_path=current_path(1,1:(lo(end)-1));
dir_path=strcat(parent_path,'/Basic_info/chrm_state/');


load Basic_info/chromosome_info.mat;
%%Path for the chromHMM data

listing = dir(dir_path);

cou=0;
for i=1:length(listing)   %%obtain all the names of chromHMM data
    lo=strfind(listing(i).name,'E');
    if ~isempty(lo)
        cou=cou+1;
        name{cou,1}=listing(i).name;
    end
end


n=length(name);


ma=zeros(15,n);
chrm_state=zeros(15,n);

for i=1:n
    te=strcat('Basic_info/chrm_state/',name{i,1});
    
    for j=1:22
        temp_data=peak_a(peak_a(:,1)==j,2:3); %%Hotspots in the current chromosome
        if ~isempty(temp_data)
            te2=strcat(te,'/chr');
            te2=strcat(te2,num2str(j));
            te2=strcat(te2,'_s.mat');
            load (te2); %%Load the chroHMM data for the current chromosome
            num=length(temp_data(:,1));
            for k=1:num
                t=loc_e(temp_data(k,1):temp_data(k,2),1);  %%The Chromatin states in the region of the current hotspot
                seed=unique(t);
                s_num=length(seed);
                if s_num==1
                    ma(seed,i)=ma(seed,i)+1;  %%If the Chromatin state has overlap with the current hotspots, then add 1.
                else
                    ma(seed(:),i)=ma(seed(:),i)+1;
                end
            end
            chrm_state(:,i)=chrm_state(:,i)+state_num;
        end
    end
end

save_name=strcat(out_res,'enri.mat');
save((save_name),'chrm_state','ma','name');  %%Save the overlap information

end
