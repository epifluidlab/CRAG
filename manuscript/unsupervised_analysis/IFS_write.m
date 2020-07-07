function IFS_write(data_name)
%%%%%%%%%%%%%%Use IFS to write the fragment information (txt file) into .mat file
%%%%%%%%%%%%%%data_name:The sample name

%%%%%%Add the path of the funtions used in this pipeline as workplace
current_path=pwd;
lo=strfind(current_path,'/');
parent_path=current_path(1,1:(lo(end)-1));
addpath(genpath(parent_path));

%%%%%%%%%Use txt_read.m to write all the fragment information chromose.
for i=1:22
    txt_read(data_name,i);
end


end
