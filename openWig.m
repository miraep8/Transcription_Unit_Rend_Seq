function reads = openWig(filename, varargin)
%OPENWIG opens the provided wig file and puts the contents into a 2xn array
% reads = OPENWIG(filename)
%
%   - filename - a string equivalent to the file path to the wig file to be
%   opened
%
%reads = OPENWIG(filename, 'header_sz', x, 'format_spec' f_str)
%
%   - filename - a string equivalent to the file path to the wig file to be
%   opened
%
%   - 'header_sz' a string which designates that the next value should be
%   assigned to the variable header_sz
%
%   - x the desired value of header_sz.  Defaults to 2.
%
%   - 'format_spec' a string which signals the user is passing a new value
%   for the format spec which will follow.
%  
%   -f_str a string which describes the format of the file.  Defaults to '%i\t%f'

    % unpack and assign changeable variables
    opts = containers.Map({'header_sz', 'format_spec'}, {2, '%i\t%f'});
    v = unpackVals(varargin, opts);
    header_sz = v(1);
    format_spec = v(2);
    
    %open the file and skip over the headers
    fileID = fopen(filename,'r');
    for i = 1:(header_sz+1)
       fgetl(fileID); 
    end
    
    %now read the rest of the file using the format specified
    reads = fscanf(fileID, format_spec, [2, Inf]);

end