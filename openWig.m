function reads = openWig(filename)

    % Used to open and read a wig file.  Returns a 2*length array of positions
    %and read counts
    %  *filename - the location of the wig file to be read
    %  *reads - the 2*length of wig array that contains position and read
    % count information. 
    
    %set up and read/clear the header:
    header_sz = 2;
    fileID = fopen(filename,'r');
    for i = 1:(header_sz+1)
       fgetl(fileID); 
    end
    
    %now read the rest of the file 
    formatSpec ='%i\t%f';
    reads = fscanf(fileID, formatSpec, [2, Inf]);

end