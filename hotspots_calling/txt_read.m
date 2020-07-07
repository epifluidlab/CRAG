function txt_read(path_name,loop_id)
%%%%%%%%%%%%Read the chromosome information
A=importdata('./Basic_info/chrome_info.txt');
chr_id=A.textdata;
chr_length=A.data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_name=strcat('./',path_name);
path_name=strcat(path_name,'/');
produce_file=strcat(path_name,'data_n');  %%%%%% fragment length

system(['mkdir ' produce_file]);

%%%%%%%%%%%%%%%%%%%%%%%%%%
path(path,path_name);
listing = dir(path_name);


cou=0;
for i=1:length(listing)   %%obtain all the names of the usefule files
    lo=strfind(listing(i).name,'ch');
    if ~isempty(lo)
        cou=cou+1;
        file_name{cou,1}=listing(i).name;
    end
end

len_threshold=1000;

for i=loop_id:loop_id
    chr=strcat('chr',num2str(i));
    name=strcat(chr,'.txt');
    lo=strcmpi(file_name,name);
    if (sum(lo)==0)
        disp('The chromosome does not exist!');
        return;
    end
    
    file_name=strcat(path_name,name);
    
    fileID = fopen(file_name);
    C = textscan(fileID,'%s %d %d %d');  %%import the bed file of the reads
    fclose(fileID);
    if isempty(C{1,1})
        disp('There is no reads in the chromosome!');
        return;  
    end
    clear read;
    B=C{1,1};
    [B,index]=sort(B);
    A=[C{1,2} C{1,3} C{1,4}];
    A=A(index,:);
    
    num=length(A);
    loop_end=num-1;
    read=zeros(num,2);
    flag_s=zeros(num,1);
    cou=0;
    for j=1:loop_end %%merge the reads that are from the same fragment
        if flag_s(j,1)==0
            if strcmpi(B{j,1},B{j+1,1})==1
                if A(j,2)>=0
                    w=j;
                else
                    w=j+1;
                end
                cou=cou+1;
                read(cou,1)=A(w,1);
                read(cou,2)=A(w,2);
                flag_s(j,1)=1;
                flag_s(j+1,1)=1;
            else
                if A(j,2)>=0
                    cou=cou+1;
                    read(cou,1)=A(j,1);
                    read(cou,2)=A(j,2);
                else
                    cou=cou+1;
                    read(cou,1)=A(j,3);
                    read(cou,2)=abs(A(j,2));
                end
                flag_s(j,1)=1;
            end
        end
    end
    if (flag_s(num,1)==0)
        if A(num,2)>=0
            cou=cou+1;
            read(cou,1)=A(num,1);
            read(cou,2)=A(num,2);
        else
            cou=cou+1;
            read(cou,1)=A(num,3);
            read(cou,2)=abs(A(num,2));
        end
    end
    
    if cou<num
        read((cou+1):num,:)=[];
    end
    
    %%remove the bad reads (fragment size more than 1000)
    read=read((read(:,2)<len_threshold) & (read(:,2)> 0),:);
    clear A;
    clear C;
    
    index=strcmpi(chr_id,chr);
    region_len=chr_length(index,1);
    read=read((read(:,1)+read(:,2))<=region_len,:);
    %%all the read are saving at 'read'
    
    
    [loc,dark_flag,read] = txt_data_prepare(read,chr,region_len);
    %%%using function txt_data_prepare to prepare the reads.
    
    dark_flag=int8(dark_flag);
    
    file_name=strcat(produce_file,'/');
    file_name=strcat(file_name,chr);
    file_name=strcat(file_name,'.mat');
    %%Saving the information to file 'chrid.mat'
    %loc:vector(n,1),where n is the number of bases on the current
    %chromosome. loc(i,1) is the IFS in the ith base
    %%dark_flag:vector(n,1).1:the base should be used;0:the base should not be used
    %%that is, the base in within dark regions.
    %%region_len, size of the current chromosome
    %%read: all the fragments in the current chromosome
    save((file_name),'loc','dark_flag','region_len','read','-v7.3');
end

end
