%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @File Name: XYZ2Lab.m
% @Description: Conversion from XYZ to LAB 
% @Params: 
%                    XYZ = input tristimulus ( n by 3)
%                    XYZn = reference white (1 by 3)
%
% @Author: Yongmin Park
% @Date: 10/18/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Lab = XYZ2Lab( XYZ, XYZn )
    if nargin < 1
        help mfilename
    elseif size(XYZ,2) < 3
        error(['Input argument error to: ' mfilename ]);
    elseif size(XYZn,2) < 3
        error(['Input argument error to: ' mfilename ]);
    elseif size(XYZn,1) < 1
        error(['Input argument error to: ' mfilename ]);
    end
    
    [row, col] = size(XYZ);
    
    % Assume XYZ is n by 3 matrix, and XYZn is 1 by 3 matrix
    M = [0 116 0;
        500 -500 0;
        0 200 -200];
    
    % Make a 3 by n matrix from the output of f function
    fXYZ = reshape([f(XYZ(:, 1) ./ XYZn(:,1)); 
        f(XYZ(:, 2) ./ XYZn(:,2));  
        f(XYZ(:, 3) ./ XYZn(:,3))], [row, col]);
    offset = [16; 0; 0];
    
    Lab = ((M * fXYZ') - offset)';    
    
end

function out = f(x)    
    out = (x > (24/ 116)^3).*(x.^(1/3)) + (x <= (24/ 116)^3).*((841/108).*x + (16/116));
end
