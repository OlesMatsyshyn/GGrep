function A = update_progress(filename)
    fileID = fopen(filename,'r');
    formatSpec = '%d';
    A = fscanf(fileID,formatSpec);
    fclose(fileID);
    fileID = fopen( filename, 'wt' );
    fprintf( fileID,'%d\n',A+1);
    fclose(fileID);
end