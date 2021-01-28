
% read in the data
reads_3f = openWig('C:/Users/mirae/Documents/Research/Data/wig Files/Jean_3f.wig');
reads_3r = openWig('C:/Users/mirae/Documents/Research/Data/wig Files/Jean_3r.wig');
reads_5f = openWig('C:/Users/mirae/Documents/Research/Data/wig Files/Jean_5f.wig');
reads_5r = openWig('C:/Users/mirae/Documents/Research/Data/wig Files/Jean_5r.wig');

%for testing purposes:
reads_3f = reads_3f(1:2,1:1000);
reads_3r = reads_3r(1:2,1:1000);
reads_5f = reads_5f(1:2,1:1000);
reads_5r = reads_5r(1:2,1:1000);
%---------------------------------

% set variables which determine the fit
d = 50; % assumed to be minimum distance between peaks
s = 500; % ' step size' how much to shift the window over by each time
w = 1000; %width of update of values - determines how 'local'

sum_f = addPeaks(reads_3f, reads_5f);
sum_r = addPeaks(reads_3r, reads_5r);
peaks_f = hmmFit(reads_3f, reads_5f, sum_f, d, s, w);
peaks_r = hmmFit(reads_3r, reads_5r, sum_r, d, s, w);

peaks_f


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

function sum_array = addPeaks(reads3, reads5)

    % Used to sum the reads3 and the reads5 arrays.  Will sum them
    % based on the indexes in the first column - adding those values 
    % with the same index to generate the value in the second column. 
    

    added = length(unique([reads3(1,:), reads5(1,:)]));
    sum_array = zeros(2, added);
    sum_array(1,:) = unique([reads3(1,:), reads5(1,:)]);
    count_3 = 1;
    count_5 = 1;
    for i = 1:added
        reads_3_value = 0;
        while reads3(1, count_3) <= i
            if reads3(1, count_3) == i
                reads_3_value = reads3(2, count_3);
            end
            count_3 = count_3 + 1;
        end
        reads_5_value = 0;
        while reads5(1, count_5) <= i
            if reads5(1, count_5) == i
                reads_5_value = reads5(2, count_5);
            end
            count_5 = count_5 + 1;
        end
        sum_array(2, i) = reads_3_value + reads_5_value;
    end
end

function peaks = hmmFit(reads3, reads5, c_reads, d, s, w)
    
  num_iter = floor((length(c_reads) - (d+w))/s);
  peaks = zeros(2, length(c_reads));
  peaks(1,:) = c_reads(1,:);
  for i = 1:num_iter
      b = i*s;
      e = i*s + d + w;
      r_3 = reads3(1:2,find(reads3(1,:)<b,1,'last')+1:find(reads3(1,:)>e,1)-1);
      r_5 = reads5(1:2,find(reads5(1,:)<b,1,'last')+1:find(reads5(1,:)>e,1)-1);
      c_r = c_reads(1:2,find(c_reads(1,:)<b,1,'last')+1:find(c_reads(1,:)>e,1)-1);
      peaks(2,b:e) = peaks(2,b:e) + hmmHelper(r_3, r_5, c_r,d);
  end
end

function [new_path, new_cost] =  makePath(val, add_cost, old_path)
    new_path = HmmPath;
    new_path.Cost = old_path.Cost + add_cost;
    new_path.Path = [old_path.Path, val];
    new_cost = old_path.Cost + add_cost;
end

function [s_c, e_c, i_c] = calcCost(path, sum_reads, reads3, reads5, d)
    b_ind = max(1, length(path) - d);
    cv = length(path) + 1;
    while sum_reads(1,b_ind) < sum_reads(1,cv) - d  %to account for gaps
        b_ind = b_ind + 1;
    end
    window_prior = path(b_ind:length(path));
    max_ind = find(abs(window_prior), 1, 'last');
    next_max = find(abs(window_prior(1:max_ind)), 1, 'last');
    mu =  mean(sum_reads(2,max_ind:max_ind + d));
    sig = std(sum_reads(2,max_ind:max_ind + d));
    i_c = 1/normpdf(sum_reads(2,cv), mu, sig);
    s_c = 1000000;
    e_c = 1000000;
   if ((max_ind - next_max) > d) or ((cv - max_ind) > d - (max_ind - next_max))
        z_score = abs((sum_reads(2,cv) -mu)/sig);
        if (z_score > 2) and (reads3(2,cv) > reads5(2,cv))
            e_c = 1/z_score;
        elseif (z_score > 2) and (reads5(2,cv) > reads3(2,cv))
            s_c = 1/z_score;
        end
    end
    
end


function peak_fits = hmmHelper(reads3, reads5, sum_reads, d)
    
    % This function actually performs the HMM fitting.  It tries to fit the
    % most likely 'fit' based on the underslying data.  A fit is defined as
    % the assignment of S(1), E(-1), I(0) to each place in sum_reads. 
    
    first = HmmPath;
    first.Cost = 1;
    first.Path = [];
    fitPaths = [first];
    costArray = [1];
    
    while true
        [M, ind] = min(costArray);
        next = fitPaths(ind);
        if next.pathLen() == length(sum_reads) - d
            break;
        end
        fitPaths(ind) = [];
        costArray(ind) = [];
        [s_c, e_c, i_c] = calcCost(next.Path, sum_reads, reads3, reads5, d);
        [start, start_c]  = makePath(1, s_c, next);
        [fin, fin_c]  = makePath(-1, e_c, next);
        [mid, mid_c]  = makePath(0,i_c, next);
        fitPaths = [fitPaths, start, mid, fin];
        costArray = [costArray, start_c, mid_c, fin_c];
        
    end
    peak_fits = next.Path;
end


