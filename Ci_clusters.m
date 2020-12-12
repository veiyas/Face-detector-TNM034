function [Cb_cluster, Cr_cluster] = Ci_clusters(Y)
Kl = 125;
Kh = 188;
WCb = 46.97;
WHCb = 14;
WLCb = 23;
WCr = 38.76;
WHCr = 10;
WLCr = 20;
yMin = 16;
yMax = 235;

Kl = 125;
Kh = 188;
yMin = 16;
yMax = 235;

CbClusterLow = WLCb + ((Y - yMin)*(WCb - WLCb)) / (Kl - yMin);
CbClusterHigh = WHCb + ((yMax - Y)*(WCb - WHCb)) / (yMax - Kh);

CrClusterLow = WLCr + ((Y - yMin)*(WCr - WLCr)) / (Kl - yMin);
CrClusterHigh = WHCr + ((yMax - Y)*(WCr - WHCr)) / (yMax - Kh);

Cb_cluster = CbClusterLow .* uint8(Y < Kl) + CbClusterHigh .* uint8(Kh < Y);
Cr_cluster = CrClusterLow .* uint8(Y < Kl) + CrClusterHigh .* uint8(Kh < Y);

end

