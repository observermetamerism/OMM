function [m] = cost_kmedoids(P, Q)    
    %m = seuclidean(P(:, 1), Q(:, 1)) + seuclidean(P(:, 2), Q(:, 2)) + seuclidean(P(:, 3), Q(:, 3));    
    m = seuclidean(P, Q);    
end

