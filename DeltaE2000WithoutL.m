%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @File Name: DeltaE2000.m
% @Description: Compute color difference using DeltaE2000
% @Params: 
%                    Lab1 = reference value (n by 3)
%                    Lab2 = test value (n by 3)
%                
% @Author: Yongmin Park (reference - https://en.wikipedia.org/wiki/Color_difference)
% @Date: 10/18/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function   DE00=DeltaE2000WithoutL(Lab1, Lab2)
     if nargin < 2
        help mfilename
    elseif size(Lab1,2) < 3
        error(['Input argument error to: ' mfilename ]);
    elseif size(Lab2,2) < 3
        error(['Input argument error to: ' mfilename ]);
    elseif size(Lab2,1) ~= size(Lab1,1) 
        error(['Input argument error to: ' mfilename ]);
     end
     
     LCh1 = Lab2LCh(Lab1);
     LCh2 = Lab2LCh(Lab2);
     
     deltaLp = 0; 
     lbar = (LCh1(:,1) + LCh2(:,1))./2; % Lbar = (L1 + L2) / 2;
     cbar = (LCh1(:,2) + LCh2(:,2))./2; % Cbar = (C1 + C2) / 2;
     
     ap1 = Lab1(:,2) + (Lab1(:,2) ./ 2) .* (1 - sqrt(cbar.^7 ./ (cbar.^7 + 25^7))); % a1
     ap2 = Lab2(:,2) + (Lab2(:,2) ./ 2) .* (1 - sqrt(cbar.^7 ./ (cbar.^7 + 25^7))); % a2
     
     c1p = sqrt(ap1.^2 + Lab1(:,3).^2);
     c2p = sqrt(ap2.^2 + Lab2(:,3).^2);
     cbarp = (c1p + c2p) ./ 2;
     deltacp = c2p - c1p;
     
     h1p = H(ap1, Lab1(:,3));
     h2p = H(ap2, Lab2(:,3));
     
     deltahp = (abs(h1p - h2p) <= 180).* (h2p - h1p) + ((abs(h1p - h2p) > 180) & (h2p <= h1p)).* ((h2p - h1p) + 360) + ((abs(h1p - h2p) > 180) & (h2p > h1p)).* ((h2p - h1p) - 360);
     deltaHp = 2 .* sqrt(c1p .* c2p) .* sind(deltahp ./ 2);
     
     Hbarp = (abs(h1p - h2p) <= 180).* ((h1p + h2p) ./2) + ((abs(h1p - h2p) > 180) & (h1p + h2p < 360)).* ((h1p + h2p + 360) ./2) + ((abs(h1p - h2p) > 180) & (h1p + h2p >= 360)).* ((h1p + h2p - 360) ./2); 
     Hbarp = ((c1p == 0) | (c2p == 0)).* ((h1p + h2p) ./2) + ((c1p ~= 0) & (c2p ~= 0)).* (Hbarp + 0);
     
     T = 1 - 0.17 .* cosd(Hbarp - 30) + 0.24 .* cosd(2 .* Hbarp) + 0.32 .* cosd(3 .* Hbarp + 6) - 0.20 .* cosd(4 .* Hbarp - 63);
     
     SL = 1 + (0.015 .* (lbar - 50).^2) ./ sqrt(20 + (lbar - 50).^2);
     SC = 1 + 0.045 .* cbarp;
     SH = 1 + 0.015 .* cbarp .* T;
     RT = -2 .* sqrt((cbarp .^ 7) ./ (cbarp .^ 7 + 25 .^ 7)) .* sind(60 .* exp(-(((Hbarp - 275) ./ 25).^2)));
     
     KL = 1;
     KC = 1;
     KH = 1;
     
     DE00 = sqrt((deltaLp ./ (KL .* SL)).^2 + (deltacp ./ (KC .* SC)).^2 + (deltaHp ./ (KH .* SH)).^2 + (RT .* deltacp .* deltaHp ./ (KC .* SC) ./ (KH .* SH)));
     
end

function out = H(a, b)    
    out = atan2d(b, a);    
    out = (out >= 0).*out + (out < 0).*(out + 360);
end