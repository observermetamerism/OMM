% Calculations to duplicate the CMSs from Mark's spreadsheet
% Using Matlab 7.3.1.1266 (R2006b) on my Macbook results agree to 10^(-15)
% read in the raw data. Note first two points of docul are interpolated.
% modified by YutaAsano (2012.10.12) (change load func to importdata)

function [ LMS, trans_lens, trans_macula, sens_photopig ] = cie2006cmfsEx( age,fs, ...
    var_od_lens, var_od_macula, var_od_L, var_od_M, var_od_S, var_shft_L, var_shft_M, var_shft_S, files )

% Load files
rmd = files.rmd; 
LMSa = files.LMSa; 
docul = files.docul; 

wl = [390:5:780]';

% field size corrected macular density
pkOd_Macula = 0.485*exp(-fs./6.132) * (1+var_od_macula/100); % varied peak optical density of macula
corrected_rmd =rmd.*pkOd_Macula;

% age corrected lens/ocular media density 
if (age<=60)
    correct_lomd = docul(:,1) .* (1+0.02*(age-32)) + docul(:,2);
else
    correct_lomd = docul(:,1) .* (1.56+0.0667*(age-60)) + docul(:,2);
end
correct_lomd = correct_lomd .* (1+var_od_lens/100); % varied overall optical density of lens

% Peak Wavelength Shift
wl_shifted(:,1) = wl + var_shft_L; 
wl_shifted(:,2) = wl + var_shft_M; 
wl_shifted(:,3) = wl + var_shft_S; 

LMSa(wl>=620,3) = NaN; % avoid interpolating above 620nm
LMSa_shft(:,1) = interp1(wl_shifted(:,1),LMSa(:,1),wl,'spline'); 
LMSa_shft(:,2) = interp1(wl_shifted(:,2),LMSa(:,2),wl,'spline'); 
LMSa_shft(:,3) = interp1(wl_shifted(:,3),LMSa(:,3),wl,'spline'); 
LMSa_shft(wl>=620,3) = NaN; 

% Surpress Warning: Columns of data containing NaN values have been ignored during interpolation. 
[~, MSGID] = lastwarn(); 
if ~isempty(MSGID); warning('off', MSGID); end


% corrected LMS (no age correction)
pkOd_L = (0.38+0.54*exp(-fs/1.333)) * (1+var_od_L/100); % varied peak optical density of L-cone
pkOd_M = (0.38+0.54*exp(-fs/1.333)) * (1+var_od_M/100); % varied peak optical density of M-cone
pkOd_S = (0.30+0.45*exp(-fs/1.333)) * (1+var_od_S/100); % varied peak optical density of S-cone

alpha_lms = 0. * LMSa_shft;
alpha_lms(:,1) = 1-10.^(-pkOd_L*(10.^LMSa_shft(:,1)));
alpha_lms(:,2) = 1-10.^(-pkOd_M*(10.^LMSa_shft(:,2)));
alpha_lms(:,3) = 1-10.^(-pkOd_S*(10.^LMSa_shft(:,3)));


% this fix is required because the above math fails for alpha_lms(3,:)==0
alpha_lms(wl>=620,3) = 0; 

% Corrected to Corneal Incidence
lms_barq = alpha_lms .* repmat(10.^(-corrected_rmd-correct_lomd),1,3);

% Corrected to Energy Terms
lms_bar = lms_barq .* repmat(wl,1,3);

% normalized
% LMS = lms_bar ./ repmat(max(lms_bar),79,1); % normalize by max
LMS = lms_bar * diag( 1./sum(lms_bar).*100 , 0 ); % normalize by equal-area


% Output extra
trans_lens = 10.^(-correct_lomd); 
trans_macula = 10.^(-corrected_rmd); 
sens_photopig = alpha_lms .* repmat(wl,1,3);

end