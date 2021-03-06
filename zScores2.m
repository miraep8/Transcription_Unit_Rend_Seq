function z_score = zScores2(reads, varargin)
% ZSCORES transforms a data file into zscores based on local data
%  z_score = ZSCORES(reads)
%
%   - reads will be the raw data (in the form of a 2Xn matrix) which 
%     contains locations in the first column and raw data values in the 
%     second.
%
%   - z_score will be the thing


%------------------------------------------------------------
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
    
    l_enough_val = l_low < i - 1 && l_low < l_high;
    r_enough_val = r_high > i + 1 && r_high > r_low;
    l_poiss = (reads(1,i - step) <= val - step - int32(.1*num_include)) && mean(reads(2, l_low:l_high)) <= 10;
    r_poiss = (reads(1,i + step) >= val + step + int32(.1*num_include)) && mean(reads(2, r_low:r_high)) <= 10;
    
    if l_enough_val && std(reads(2, l_low:l_high)) ~= 0 && ~l_poiss
        l_mean = mean(reads(2, l_low:l_high));
        l_std = std(reads(2, l_low:l_high));
        l_score = (reads(2, i) - l_mean)/l_std;
    elseif l_poiss
        poiss_val = poissFind(i, val - step, val - gap, reads);
        lam = mean(poiss_val);
        lam = max(1,lam);
        l_p_z = poisscdf(reads(2,i), lam);
        l_score = icdf('Normal', l_p_z, 0, 1);
    end
    if  r_enough_val && std(reads(2, r_low:r_high)) ~= 0
        r_mean = mean(reads(2, r_low:r_high));
        r_std = std(reads(2, r_low:r_high));
        r_score = (reads(2, i) - r_mean)/r_std;
    elseif r_poiss
        poiss_val = poissFind(i, val + gap, val + step, reads);
        lam = mean(poiss_val);
        lam = max(1,lam);
        r_p_z = poisscdf(reads(2,i), lam);
        r_score = icdf('Normal', r_p_z, 0, 1);
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

function vals = poissFind(cur_ind, low_v, high_v, reads)

vals = zeros(1,high_v - low_v);
for i = low_v:high_v
    while reads(1,cur_ind) < i
        cur_ind = cur_ind + 1;
    end
    if reads(1,cur_ind) == i
        vals(i) = reads(2, cur_ind);
        cur_ind = cur_ind + 1;
    end
end

end

function index = adjustDown(cur_val, target_val, reads)
    while reads(1,cur_val) > target_val
        cur_val = cur_val - 1;        
    end
    index = cur_val;
end

function index = adjustUp(cur_val, target_val, reads)
    while reads(1,cur_val) < target_val
        cur_val = cur_val + 1;        
    end
    index = cur_val;
end




