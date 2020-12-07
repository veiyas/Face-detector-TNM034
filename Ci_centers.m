function [Cb_center, Cr_center] = Ci_centers(Y)
Kl = 125;
Kh = 188;
yMin = 16;
yMax = 235;

if Y < Kl
    Cb_center = 108 + ((Kl-Y)*10)/(Kl - yMin);
    Cr_center = 154 - ((Kl-Y)*10)/(Kl - yMin);
elseif Kh <= Y
    Cb_center = 108 + ((Y-Kh)*10)/(yMax - Kh);
    Cr_center = 154 + ((Y-Kh)*22)/(yMax - Kh);
end


end

