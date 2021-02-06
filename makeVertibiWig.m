function wig_track = makeVertibiWig(z_scores, varargin)
%MAKEVERTIBIWIG take a 2xn array of z scores and applies the Vertibi
%algorithm to assign hidden states to the observed data (z_scores)
%
%wig_track = MAKEVERTIBIWIG(z_scores)
%  
%  - z_scores a 2xn array (position, z_transformed data) of z transformed data to be fit
%
%  - wig_track - a 2xn array (position, hidden state assignment) of the
%  vertibi algorithms fit of the raw data (using the HMM)
%
%  - optional variables (to use an optional variable - pass in the name as
%  a string, followed by your desired value i.e MAKEVERTIBIWIG(z_scores,
%  'p' .1))
%
%          - 'p' the probability of transferring to a the peak state while
%          at the internal state.  Should be equal to 1/[length of
%          transcription unit] (assuming a geometric distribution). Default
%          is 1/2000. 
%
%          - 'cutoff' the value of the z score below which the probability
%          that a peak generated that z-score is zero. Default is 7. 

    opts = containers.Map({'p', 'cutoff'}, {1/2000, 7});
    v = unpackVals(varargin, opts);
    p = v(1);
    cutoff = v(2);
    
    states = [1,100]; % I, P  - the states this is what wig will write
    trans_M = [[(1-p),(p)];[1,0]]; %transition probability matrix
    wig_track = zeros(2, length(z_scores));
    wig_track(1,:) = z_scores(1,:);
    
    T1 = zeros(length(states), length(z_scores));
    T2 = zeros(length(states), length(z_scores));
    T1(:,1) = ones(length(states),1);
    norm_val = normcdf(cutoff, cutoff+5, 10); %to normalize peak emission probs

    %Vertibi Algorithm:
    for i = 2:length(z_scores)
        p_I = normpdf(z_scores(2,i));
        p_P = 0;
        if z_scores(2,i) > cutoff
            p_P = normpdf(z_scores(2,i),5,10)/(1-norm_val);
        end
        probs = [p_I, p_P];  %emission probabilities
            
        for j = 1:length(states)
            paths = zeros(length(states), 1);
            for k = 1:length(states)
                if (T1(k, i-1) == -Inf) || (trans_M(k,j) == 0) || (probs(j) == 0)
                    paths(k) = -Inf;
                else
                    paths(k) = T1(k, i-1) + log(trans_M(k, j)) + log(probs(j));  
                end
            end  
            [~, T2(j,i)] = max(paths);  
            T1(j,i) = paths(T2(j,i));
        end
    end
    
    max_inds = zeros(length(wig_track), 1);
    [~, max_inds(length(wig_track))] = max((T1(:, length(T1))));
    
    wig_track(2,length(wig_track)) = states(max_inds(length(wig_track)));
    for m = flip(2:length(wig_track))
        max_inds(m-1) = T2(max_inds(m), m);
        wig_track(2,m-1) = states(max_inds(m-1));
    end 
    
    
end