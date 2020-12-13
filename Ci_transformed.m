function [Cb_xform, Cr_xform] = Ci_transformed(Cb, Cr, CbCluster, CrCluster, CbCenter, CrCenter)
%Ci_transformed makes the final transform to luma-independant YCbCr colors
%   Detailed explanation goes here
Wcb = 46.97;
Wcr = 38.76;
Kb = 184;
Kr = 168;

Cb_xform = (Cb - CbCenter).*(Wcb ./ CbCluster) + Kb;
Cr_xform = (Cr - CrCenter).*(Wcr ./ CrCluster) + Kr;
end

