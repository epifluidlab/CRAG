function IFS_write(data_name)
%%%%%%%%%%%%%%Use IFS to write the fragment information (txt file) into .mat file
%%%%%%%%%%%%%%data_name:The sample name

%%%%%%%%%Use txt_read.m to write all the fragment information chromose.
for i=1:22
    txt_read(data_name,i);
end


end