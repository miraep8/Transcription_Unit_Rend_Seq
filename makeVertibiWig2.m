function wig_track = makeVertibiWig2(z_scores)

    p = 1/2000; % p of a new peak while in i.  2000 is approx the length of t.u.
    cutoff = 7;
    states = [1,100]; % I, P  - the states
    trans_M = [[(1-p),(p)];[1,0]];
    wig_track = zeros(2, length(z_scores));
    wig_track(1,:) = z_scores(1,:);
    
    T1 = zeros(length(states), length(z_scores));
    T2 = zeros(length(states), length(z_scores));
    T1(:,1) = ones(length(states),1);
    notyet = 0;
    for i = 2:length(z_scores)
        p_I = normpdf(z_scores(2,i));
        p_P = 0;
        if z_scores(2,i) > cutoff
            norm_val = normcdf(cutoff, cutoff+5, 10);
            p_P = normpdf(z_scores(2,i),5,10)/(1-norm_val);
        end
        probs = [p_I, p_P];
            
        for j = 1:length(states)
            paths = zeros(length(states), 1);
            for k = 1:length(states)
                if (T1(k, i-1) == -Inf) || (trans_M(k,j) == 0) || (probs(j) == 0)
                    paths(k) = -Inf;
                else
                    paths(k) = T1(k, i-1) + log(trans_M(k, j)) + log(probs(j));  
                end
            end  
            if (notyet == 0) && (j ==1) && (sum(paths == -Inf) > 1)
                i
                notyet = 1;
                p = probs
                tM = trans_M
                T1(:,i-10:i)
            end
            [~, T2(j,i)] = max(paths);  
            T1(j,i) = paths(T2(j,i));
        end
    end
    
    max_inds = zeros(length(wig_track), 1);
    [n, max_inds(length(wig_track))] = max((T1(:, length(T1))));
    n
    wig_track(2,length(wig_track)) = states(max_inds(length(wig_track)));
    for m = flip(2:length(wig_track))
        max_inds(m-1) = T2(max_inds(m), m);
        wig_track(2,m-1) = states(max_inds(m-1));
    end 
    
    
end