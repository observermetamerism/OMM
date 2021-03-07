function [ varParam ] = fnc_MonteCarloParam( stdDevAllParam )
%UNTITLED Summary of this function goes here

%{ 

list_paramNames = {'od_lens', 'od_macula', 'od_L', 'od_M', 'od_S', 'shft_L', 'shft_M', 'shft_S'}; 

stdDevAllParam(1) = 12.5; 
stdDevAllParam(2) = 27; 
stdDevAllParam(3) = 15; 
stdDevAllParam(4) = 15; 
stdDevAllParam(5) = 15; 
stdDevAllParam(6) = 2; 
stdDevAllParam(7) = 1.5; 
stdDevAllParam(8) = 1.25; 

%}

varParam = nan(size(stdDevAllParam)); 
for k = 1:length(stdDevAllParam)
    varParam(k) = stdDevAllParam(k) * randn(1); 
    
    % limit varAllParam so that it doesn't create negative val for lens, macula, pkod_LMS
    if k == 1 || 2|| 3 || 4 || 5
        if varParam(k)<-100
            varParam(k) = -100; 
        end
    end
end



end

