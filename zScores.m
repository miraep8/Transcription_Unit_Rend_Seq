
function z_score = zScores(reads, varargin)
% ZSCORES transforms a data file into zscores based on local data
%  z_score = ZSCORES(reads)
%
%   - reads will be the raw data (in the form of a 2Xn matrix) which 
%     contains locations in the first column and raw data values in the 
%     second.
%
%   - z_score will be the z_score transformed data (ie each value will now
%   be transformed by assuming that the values around it (to the left and
%   right) constitute normal distributions. Each value is replaced with the
%   minimum of its z_score from the left and right distributions. 
%
%  - optional variables (to use an optional variable - pass in the name as
%  a string, followed by your desired value i.e z_s = ZSCORES(reads,
%  'gap' 5))
%
%          - 'gap' is the number of spaces immediately around the
%          nucleotide (in terms of position space not space in the array)
%          to ignore when calculating the normal distributions. 
%
%          - 'num_include' is the number of values below the gap to use in
%          calculating the normal distribution.  This is number of values
%          in position space - so may actually end up being fewer values if
%          reads are sparse.

opts = containers.Map({'gap', 'num_include'}, {5, 50});
v = unpackVals(varargin, opts);
gap = v(1);
num_include = v(2);
step = gap + num_include;

z_score = zeros(2, length(reads) - 2*(step));
z_score(1,:) = reads(1,(step + 1):length(reads) - (step));

for i = (step + 1):(length(reads) - (step))
    val = reads(1,i);
    l_low = adjustUp(i - step, val-step, reads);
    l_high = adjustDown(i - gap, val - gap, reads);
    r_low = adjustUp(i + gap, val + gap, reads);
    r_high = adjustDown(i + step, val + step, reads);
    if l_low < i - 1 && l_low < l_high && std(reads(2, l_low:l_high)) ~= 0
        l_mean = mean(reads(2, l_low:l_high));
        l_std = std(reads(2, l_low:l_high));
        l_score = (reads(2, i) - l_mean)/l_std;
    else
        l_score = reads(2,i);
    end
    if r_high > i + 1 && r_high > r_low && std(reads(2, r_low:r_high)) ~= 0
        r_mean = mean(reads(2, r_low:r_high));
        r_std = std(reads(2, r_low:r_high));
        r_score = (reads(2, i) - r_mean)/r_std;
    else
        r_score = reads(2,i);
    end
    
    z_score(2,i - step) = r_score;
    if min(abs(l_score), abs(r_score)) == abs(l_score)
        z_score(2,i - step) = l_score;
    end
    if abs(z_score(2,i - step)) > 15
        z_score(2,i - step) = 15*sign(z_score(2,i - step));
    end
    
end

end

function index = adjustDown(cur_ind, target_val, reads)
%helper function.  Is neccesary because there are sometimes gaps in the
%position space.  Will adjust the cur_ind down until its corresponding
%position in reads is below target_val

    while reads(1,cur_ind) > target_val
        cur_ind = cur_ind - 1;        
    end
    index = cur_ind;
end

function index = adjustUp(cur_ind, target_val, reads)
%helper function.  Is neccesary because there are sometimes gaps in the
%position space.  Will adjust the cur_ind up until its corresponding
%position in reads is above target_val

    while reads(1,cur_ind) < target_val
        cur_ind = cur_ind + 1;        
    end
    index = cur_ind;
end



