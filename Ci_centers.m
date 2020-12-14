function [Cb_center, Cr_center] = Ci_centers(Y)
Kl = 125;
Kh = 188;
yMin = 16;
yMax = 235;

CbCenterLow = 108 + (Kl - Y) * 10 / (Kl - yMin);
CbCenterHigh = 108 + (Y - Kh) * 10 / (yMax - Kh);

CrCenterLow = 154 - (Kl - Y) * 10 / (Kl - yMin);
CrCenterHigh = 154 + (Y - Kh) * 22 / (yMax - Kh);

Cb_center = CbCenterLow .* uint8(Y < Kl) + CbCenterHigh .* uint8(Kh < Y);
Cr_center = CrCenterLow .* uint8(Y < Kl) + CrCenterHigh .* uint8(Kh < Y);

% if Y < Kl
%     Cb_center = 108 + ((Kl-Y)*10)/(Kl - yMin);
%     Cr_center = 154 - ((Kl-Y)*10)/(Kl - yMin);
% elseif Kh <= Y
%     Cb_center = 108 + ((Y-Kh)*10)/(yMax - Kh);
%     Cr_center = 154 + ((Y-Kh)*22)/(yMax - Kh);
% else
%     Cb_center = 0;
%     Cr_center = 0;
% end


end

