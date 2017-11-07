function [ fid_file, fid_txt ] = binary_saver( name_file, data, num_file )

    fid_file = fopen([name_file, '_' ,sprintf('%3.3d',num_file), '.bin'], 'w');
    txt_name = [name_file, '_spec'];
    fid_txt = fopen([txt_name, '_', sprintf('%3.3d',num_file), '.txt'], 'wt');
    
    data_size = size(data);
    
    fprintf(fid_txt, '%d,%d,float32', data_size(1), data_size(2));
    fwrite(fid_file, data, 'float32');
    fclose(fid_file);
    fclose(fid_txt);

end

