function [loc,dark_flag,read]=txt_data_prepare(read,chr,region_len)

[data,txt]=xlsread('Blacklist_hg19.xlsx');
loc_black=strcmpi(txt(:,1),chr);  %%get the black region

path(path,'Basic_info/mappability/');
m_name='Basic_info/mappability/';

m_name=strcat(m_name,chr);
m_name=strcat(m_name,'_m.mat');  %%get the mappability scores
load (m_name);

if sum(loc_black)>0
    black_region=data(loc_black,1:2);
    black_region=sortrows(black_region,1);
    m=length(black_region(:,1));
    dark_flag=ones(region_len,1);
    
    for i=1:m
        dark_flag(black_region(i,1):black_region(i,2),1)=0;
    end
else
    dark_flag=ones(region_len,1);
end
dark_flag = (dark_flag & (map_b~=0));
%%%%get the dark_flag: vector(n,1).1:the base should be used;0:the base should not be used
%%that is, the base in within dark regions or with zero mappability score.


point=read(:,1)+fix(read(:,2)/2); %% median point for each fragment

%%Only the fragments don't locate within dark regions or with zero mappability score were used.
read=read(dark_flag(point)==1,:);
point=point(dark_flag(point)==1,1);

total_read=length(read);
total_len=sum(read(:,2));
loc=zeros(region_len,1);
n=length(read);
for i=1:n
    if dark_flag(point(i,1),1)==1
       %%IFS for all the base in the chrrent chromosome
        loc(point(i,1),1)=loc(point(i,1),1)+1;
        loc(point(i,1),1)=loc(point(i,1),1)+((read(i,2)./total_len)*total_read);
    end
end

end
