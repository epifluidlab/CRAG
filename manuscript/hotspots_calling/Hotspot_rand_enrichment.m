function Hotspot_rand_enrichment(data_name)
%%%%%%%Calculate the overlap number of random hotspots (The same size with true hotspots) for each chromHMM, which
%%%%%%%would laters used to chromHMM enrichment
%%%%%%%%%data_name: The file name of the data set. i.e. BH01

file_path = strcat('./',data_name);
out_res= strcat(file_path,'/result_n/');
file_path = strcat(file_path,'/result_n/peak_all.mat');


load (file_path);   %%load all the hotspots data

path(path,'./Basic_info/');


load ./Basic_info/chromosome_info.mat;

path(path,'./Basic_info/chrm_state/');   %%chr1_m.mat

listing = dir('./Basic_info/chrm_state/');


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
    te=strcat('./Basic_info/chrm_state/',name{i,1});
    for j=1:22
        r_data=peak_a(peak_a(:,1)==j,2:3); %%Hotspots in the current chromosome
        if ~isempty(r_data)
            te2=strcat(te,'/chr');
            te2=strcat(te2,num2str(j));
            te2=strcat(te2,'_s.mat');
            load (te2); %%Load the chroHMM data for the current chromosome
            p_num=length(r_data(:,1));
            for r=1:10  %%Repeat the ramdom processes 10 times.
                temp_data=zeros(p_num,2);
                
                for w=1:p_num
                    rand_l=randi(chrm_len(j,1)-(r_data(w,2)-r_data(w,1)),1,1);  %%Produce the same numerb of random hotspots, which were also with the same length with the true hotspots
                    temp_data(w,1)=rand_l;
                    temp_data(w,2)=rand_l+(r_data(w,2)-r_data(w,1));
                end
                
                num=length(temp_data(:,1));
                for k=1:num
                    t=loc_e(temp_data(k,1):temp_data(k,2),1); %%The Chromatin states in the region of the current random hotspot
                    seed=unique(t);
                    seed=setdiff(seed,0);  %%remove state 0
                    s_num=length(seed);
                    if s_num~=0
                        ma(seed(:),i)=ma(seed(:),i)+1;  %%If the Chromatin state has overlap with the current random hotspot, then add 1.
                    end
                end
            end
            chrm_state(:,i)=chrm_state(:,i)+state_num;
        end
    end
end
ma_c=ceil(ma./10);  %% Take the average overlap values in the 10 repeats as the final overlap value for each chromHMM.

save_name=strcat(out_res,'enri_randm.mat');   %%Save the overlap information for the random hotspots.
save((save_name),'chrm_state','ma_c','name');

end