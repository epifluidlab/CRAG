function myth_signal_hotspot(hotspot_file,myth_file)
%%%Calculate the DNA mythelation levels around hotspots
%%hotspot_file:the name of the hotspot sample:i.e.BH01
%%myth_file:the sample name of DNA mythelation levels, i.e.PNAS_healthy
%%By default, the myth_file should be located in './Basic_info/myth_path'
%% such as './Basic_info/myth_path/PNAS_healthy/'.

peak_file=strcat(hotspot_file,'/result_n/peak_all');
load (peak_file);

file_path =strcat('Basic_info/myth_path/',myth_file);
out_res='Basic_info/myth_path/';


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

data=zeros(201,2);
if (cou~=0)
    num=length(name);
    for ii=1:num
        te=name{ii,1};
        lo=strfind(te,'.');
        chr=str2double(te(1,4:(lo-1)));
        
        t_name=strcat(file_path,'\');
        read_file=strcat(t_name,name{ii,1});
        
        load (read_file);
        peak_small=peak_a(peak_a(:,1)==chr,2:3);
        if (~isempty(peak_small))
            n=length(peak_small(:,1));
            for i=1:n
                base=ceil((peak_small(i,1)+peak_small(i,2))/2);
                if ((base-2000)>=1) && ((base+2000) <= length(loc))
                    temp=loc((base-10):(base+10),1);
                    if (sum(temp~=-1)~=0)
                        data(101,1)=data(101,1)+mean(temp(temp~=-1,1));
                        data(101,2)=data(101,2)+1;
                    end
                    
                    for j=1:100
                        a=(base+(j-1)*20+11);
                        b=(base+j*20+10);
                        temp = loc(a:b,1);
                        if (sum(temp~=-1)~=0)
                            data(101+j,1)=data(101+j,1)+mean(temp(temp~=-1,1));
                            data(101+j,2)=data(101+j,2)+1;
                        end
                    end
                    
                    for j=1:100
                        a=(base-(j-1)*20-11);
                        b=(base-j*20-10);
                        temp= loc(b:a,1);
                        if (sum(temp~=-1)~=0)
                            data(101-j,1)=data(101-j,1)+mean(temp(temp~=-1,1));
                            data(101-j,2)=data(101-j,2)+1;
                        end
                    end
                end
            end
        end
    end
    
    t_file=strcat(hotspot_file,'_');
    t_file=strcat(t_file,myth_file);
    t_file=strcat(t_file,'_DNAmyth.mat');
    x_val(:,1)=(1:20:4001)';
    x_val(:,1)=x_val(:,1)-2001;
    y_val=data(:,1)./data(:,2);
    name_data=strcat(out_res,t_file);
    save((name_data),'x_val','y_val','-v7.3');
    %%Save the DNA mythelation level values.
end

end

