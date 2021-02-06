function writeWig(wig_track, wig_file_name)
% WRITEWIG Saves a 2xn data matrix into a wig file at a specified filename
%  WRITEWIG(wig_track, wig_file_name)
%
%   - wig_track should be a 2 by n data matrix where the first column is 
%   read location (in ascending order) and the second column is the data 
%   value at that point. 
%
%   - wig_file_name should be the file-path of your new desired file for 
%   the data in wig_track

    wig_track(:,wig_track(1,:) < 1) = [];
    fileID = fopen(wig_file_name,'w');
    fprintf(fileID, 'track type=wiggle_0\n');
    fprintf(fileID, 'variableStep chrom=NC_000964.3\n');
    fprintf(fileID,'%d\t%d\n',wig_track);
    fclose(fileID);   
    
end