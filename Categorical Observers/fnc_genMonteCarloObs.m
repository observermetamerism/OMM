function [LMS_All, var_age, vAll ] = fnc_genMonteCarloObs( n_population, list_Age, fs)

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%

list_paramNames = {'od_lens', 'od_macula', 'od_L', 'od_M', 'od_S', 'shft_L', 'shft_M', 'shft_S'}; 

stdDevAllParam(1) = 19.1; 
stdDevAllParam(2) = 37.2; 
stdDevAllParam(3) = 17.9; 
stdDevAllParam(4) = 17.9; 
stdDevAllParam(5) = 14.7; 
stdDevAllParam(6) = 4.0; 
stdDevAllParam(7) = 3.0; 
stdDevAllParam(8) = 2.5; 

% Scale down StdDev by scalars optimized using Asano's 75 observers collected in Germany
stdDevAllParam(1:2) = stdDevAllParam(1:2)*0.98; 
stdDevAllParam(3:end) = stdDevAllParam(3:end)*0.50; 

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%


% Get Normally-distributed Physiological Factors
vAll = nan(n_population, 8); 
for k = 1:n_population
    [ vAll(k,:) ] = fnc_MonteCarloParam( stdDevAllParam ); 
end
% stat(vAll)

% Generate Random Ages with the same probability density distribution as color matching experiment
sz_interval = 1; 
list_AgeRound = round( list_Age./sz_interval ) .* sz_interval; 
var_age = randsample(unique(list_AgeRound), n_population, true, hist(list_AgeRound,unique(list_AgeRound))); 

% 
vAll(1, :) = zeros(1, 8);
var_age(1) = 38;

files.rmd = importdata('cie2006_RelativeMacularDensity.txt');
files.LMSa = importdata('cie2006_Alms.txt');
files.docul = importdata('cie2006_docul.txt');

LMS_All = nan(79,3,n_population); 

for k = 1:n_population
[ t_LMS, ~, ~, ~ ] = cie2006cmfsEx( var_age(k),fs, ...
    vAll(k,1), vAll(k,2), vAll(k,3), vAll(k,4), vAll(k,5), vAll(k,6), vAll(k,7), vAll(k,8), files ); 

LMS_All(:,:,k) = t_LMS; 
% trans_lens_All(:,k) = trans_lens; 
% trans_macula_All(:,k) = trans_macula; 
% sens_photopig_All(:,:,k) = sens_photopig; 

end


end

