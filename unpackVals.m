function vals = unpackVals(var_input, opts)
% UNPACKVALS will take the variable arguments passed to varargin and unpack
% them according to the string values in opts, checking that they are the
% correct type.  Varargin values must be passed as a string value pair for
% this to work.
%   vals = UNPACKVALS(varargin, opts)
%
%       - var_input is a cell array of the variable user inputs.  It should
%       be of the form {string}, {value}, {string}, {value}.....
%
%       - opts is a map from characters to values.  The values should be
%       the default values for all optional input variables. The keys will
%       be the names of these variables as strings.  When users call the
%       original function they should enter optional variables in the form
%       'var_name', var_val) when calling the original function - after all
%       required variables. 
%
%       -vals will be an array of unpacked values.  If a user passed a
%       value for one of the 

vals = zeros(1, length(opts));
count = 1;
for v = keys(opts)
    find_val = strcmp(var_input, v);
    if any(find_val)
        new_val = var_input{find(find_val) + 1};
        if isa(new_val, class(opts(char(v))))
            vals(count) = new_val;
        else
            disp('Incorrect variable type passed!')
            disp('Variable ' + opts(v) + ' should be of type ' + str(class(defaults(v))))
            disp('Instead variable of type ' + str(class(new_val)) + ' was passed')
        end
    else
        vals(count) = opts(char(v));
    end
    count = count + 1;
end