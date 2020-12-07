function [WCb,WCr] = Ci_weights(Y)
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

if Y < Kl
    WCb = WLCb + ((Y - yMin)*(WCb - WLCb)) / (Kl - yMin);
    WCr = WLCr + ((Y - yMin)*(WCr - WLCr)) / (Kl - yMin);
elseif Kh < Y
    WCb = WHCb + ((yMax - Y)*(WCb - WHCb)) / (yMax - Kh);
    WCr = WHCr + ((yMax - Y)*(WCr - WHCr)) / (yMax - Kh);
end

end

