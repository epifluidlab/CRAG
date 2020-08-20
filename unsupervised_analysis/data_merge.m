function data_merge(in_path_name,out_path,file_list,loop_id)
%%%%%%%%%%%%Merge the IFS score of multi-samples and output the merged IFS
%%%%%%%%%%%%score
A=importdata('Basic_info/chrome_info.txt');
chr_id=A.textdata;
chr_length=A.data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_name=out_path;
path_name=strcat(path_name,'/');

produce_file=strcat(path_name,'data_n');  %%%%%% fragment length
sub_path='/data_n/';

if ~exist(produce_file,'dir')
system(['mkdir ' produce_file]);
end

data_path=in_path_name;
data_path=strcat(data_path,'/');


for i=loop_id:loop_id
    
    chr=strcat('chr',num2str(i));
    id_loc=strcmpi(chr_id(:,1),chr);
    loc_a=zeros(chr_length(id_loc,1),1);
    for j=1:length(file_list)
        file_name=strcat(data_path,file_list{j,1});
        file_name=strcat(file_name,sub_path);
        file_name=strcat(file_name,chr);
        file_name=strcat(file_name,'.mat');
        if exist(file_name,'file')==0
            er_info=strcat(file_name,'does not exist!');
            disp(er_info);
            return;
        end
        load (file_name);
        loc_a=loc_a+loc;
    end
    clear loc;
    loc=loc_a;
    file_name=strcat(produce_file,'/');
    file_name=strcat(file_name,chr);
    file_name=strcat(file_name,'.mat');
    save((file_name),'loc','dark_flag','region_len','-v7.3');
end

end
