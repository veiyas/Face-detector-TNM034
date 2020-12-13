function [NLYCbCr] = nlYCbCr(YCbCr)
Y = YCbCr(:,:,1);
Cb = YCbCr(:,:,2);
Cr = YCbCr(:,:,3);

% Constants attained from Face Detection in Color Images
Kl = 125;
Kh = 188; 

% Matrixes representing low/high
Yl = uint8(Y < Kl);
Yh = uint8(Y > Kh);
Yi = uint8(Y >= Kl) .* uint8(Y <= Kh);

[CbCenter, CrCenter] = Ci_centers(Y);
[CbCluster, CrCluster] = Ci_clusters(Y);
[CbXformed, CrXformed] = Ci_transformed(Cb, Cr, CbCluster, CrCluster, CbCenter, CrCenter);

% combine the results of the different threshold images
NLYCbCr = Y;
NLYCbCr(:, :, 2) = CbXformed.* uint8((Yl+Yh) >= 1) + Cb .* Yi;
NLYCbCr(:, :, 3) = CrXformed.* uint8((Yl+Yh) >= 1) + Cr .* Yi;
end