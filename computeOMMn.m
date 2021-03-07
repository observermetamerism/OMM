%%%%%%%%%%%%%%%%%%%%%%%%%%
% computeOMMn
% ** params
% refSPD: The SPDs of a reference display
% testSPD: The SPDs of a test display
% ind_xyz_CMFs: The xyz-like CMFs of a given standard observer
% color: The color for evaluation
function [omm] = computeOMMn(refSPD, testSPD, ref_xyz_CMFs, color)

     XYZref = 683 .* (ref_xyz_CMFs' * sum(refSPD)');
     XYZtest = 683 .* (ref_xyz_CMFs' * sum(testSPD)');     
     omm = DeltaE2000WithoutL(XYZ2Lab(XYZref', color), XYZ2Lab(XYZtest', color));
    
end
