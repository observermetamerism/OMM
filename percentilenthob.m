function [perc] = percentilenthob(M, n)
    [row, col] = size(M);
    
    threshold = int32(n .* row);   
    
    perc = zeros(col, 1);
    
    for i = 1:col
        eo = sort(M(:, i));
        perc(i) = eo(threshold);
    end  
    
end

