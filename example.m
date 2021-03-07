clc;
clear;
close all;

stdOb = load('CIE1931_Standard_Observer.mat'); % CIE 1931 2-deg standard observer
indObs = load('./CMFs/Cat_2deg_CMFs_N=100.mat'); % 100 Categorical observers

% CIE D65 neutral pattern
neutral = [95.0455927000000;100;108.905775100000] .* 0.5;
wavelength = 390:1:780;

[~, ~, num_of_indObs] = size(indObs.xyz_CMFs);
indObs.xyz_CMFs_interp = zeros(length(wavelength), 3, num_of_indObs);

for x = 1:num_of_indObs
    iobs = squeeze(indObs.xyz_CMFs(:, :, x));
    iobs_interp = zeros(length(wavelength), 3);
    
    indObs.xyz_CMFs_interp(:, :, x) = interp1(390:5:780,  iobs, wavelength, 'pchip');
end

refDisplay = load('./Example displays/01_Rec709.mat');
testDisplay = load('./Example displays/08_Rec2020_100.mat');

% Look at Equation (6). 
M_ref_std = (683 .* stdOb.xyz_CMFs') * refDisplay.SPD';
RGB_ref_std = M_ref_std \ neutral;
SPD_ref_std = refDisplay.SPD .* RGB_ref_std;

OMMn = zeros(num_of_indObs, 1);

for x = 1:num_of_indObs    
    % Look at Equation (7) ~ (9)
    ind_xyz_CMFs = squeeze(indObs.xyz_CMFs_interp(:, :, x));
    M_test_ind = (683 .* ind_xyz_CMFs') * testDisplay.SPD';
    XYZ_ref_ind = (683 .* ind_xyz_CMFs') * sum(SPD_ref_std)';
     
    RGB_test_ind = M_test_ind \ XYZ_ref_ind;
    SPD_test_ind = testDisplay.SPD .* RGB_test_ind;      
    
    % Look at Equation (10)
    OMMn(x) = computeOMMn(SPD_ref_std, SPD_test_ind, stdOb.xyz_CMFs, neutral');    
end

% Look at Equation (12), 90th percentile observer
perc90 = percentilenthob(OMMn, 0.90)