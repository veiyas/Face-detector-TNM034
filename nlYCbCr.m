function [NLYCbCr] = nlYCbCr(YCbCr)
    Y = YCbCr(:, :, 1);
    Cb = YCbCr(:, :, 2);
    Cr = YCbCr(:, :, 3);
    
    % neat constants from appendix A of "Face Detection in Color Images"
    % spread values
    Wcb = 46.97;
    WLcb = 23;
    WHcb = 14;
    Wcr = 38.76;
    WLcr = 20;
    WHcr = 10;
    
    Kl = 125;
    Kh = 188; 
    Ymin = 16;
    Ymax = 235;
    
    % binary threshold images
    Yl = uint8(Y < Kl);
    Yh = uint8(Y > Kh);
    Yi = uint8(Y >= Kl) .* uint8(Y <= Kh);
    % calc center values Cb, Equation 7
    Ba = 108;
    Bb = 118;
    CbCenterl = Ba + (Kl - Y) * (Bb - Ba) / (Kl - Ymin);
    CbCenterh = Ba + (Y - Kh) * (Bb - Ba) / (Ymax - Kh);
   
    CbCenter = CbCenterl .* Yl + CbCenterh .* Yh;

    % calc center values Cr, Equation 8
    Ra = 154;
    Rb = 144;
    Rc = 132;
    CrCenterl = Ra - (Kl - Y) * (Ra - Rb) / (Kl - Ymin);
    CrCenterh = Ra + (Y - Kh) * (Ra - Rc) / (Ymax - Kh);
    
    CrCenter = CrCenterl .* Yl + CrCenterh .* Yh;

    % calculate spread of b equation 6
    SpreadBl = clusterSpreadL(WLcb, Y, Ymin, Wcb, Kl);
    SpreadBh = clusterSpreadH(WHcb, Y, Ymax, Wcb, Kh);

    SpreadB = SpreadBl .* Yl + SpreadBh .* Yh;

    % calculate spread of r equation 6
    SpreadRl = clusterSpreadL(WLcr, Y, Ymin, Wcr, Kl);
    SpreadRh = clusterSpreadH(WHcr, Y, Ymax, Wcr, Kh);

    SpreadR = SpreadRl .* Yl + SpreadRh .* Yh;
    
    % actually calculate the new Cb and Cr, equation 5
    krConst = 168;
    kbConst = 184;
    CPrimB = cPrim(Cb, Wcb, SpreadB, CbCenter, kbConst);
    CPrimR = cPrim(Cr, Wcr, SpreadR, CrCenter, krConst);
    % combine the results of the different threshold images
    CPrimBi = Cb .* Yi;
    CprimB = CPrimB.*(Yl+Yh) + CPrimBi;
    
    CPrimRi = Cr .* Yi;
    CprimR = CPrimR.*(Yl+Yh) + CPrimRi;
    % merge into one image
    NLYCbCr = Y;
    NLYCbCr(:, :, 2) = CprimB;
    NLYCbCr(:, :, 3) = CprimR;
end

% Equation 5
function [res] = cPrim(C, Wc, clusterSpread, Ccenter, K)
    res = C - Ccenter;
    res = res * Wc ./ clusterSpread;
    res = res + K;
end

% Equation 6
function [res] = clusterSpreadL(WLc, Y, Ymin, Wc, Kl)
    res = WLc + ((Y - Ymin) * (Wc - WLc)) / (Kl - Ymin);
end

function [res] = clusterSpreadH(WHc, Y, Ymax, Wc, Kh)
    res = WHc + ((Ymax - Y) * (Wc - WHc)) / (Ymax - Kh);
end






















% imSize = size(YCbCr);
% Kl = 125;
% Kh = 188;
% WCb = 46.97;
% WCr = 38.76;
% 
% NLYCbCr = zeros(imSize(1), imSize(2), imSize(3));
% NLYCbCr = uint8(NLYCbCr);
% NLYCbCr(:,:,1) = YCbCr(:,:,1);
% for i = 1:imSize(1)
%     for j = 1:imSize(2)
%         Y = YCbCr(i,j,1);
%         
%         if Y < Kl || Kh < Y
%             [CbC, CrC] = Ci_centers(Y); % Centers
%             [CbCK, CrCK] = Ci_centers(Kh); % Centers with Kh as input
%             [CbW, CrW] = Ci_weights(Y); % Weights
%             NLCb = (Y - CbC)*(WCb / CbW) + CbCK;
%             NLCr = (Y - CrC)*(WCr / CrW) + CrCK;
%             NLYCbCr(i,j,2) = uint8(NLCb);
%             NLYCbCr(i,j,3) = uint8(NLCr);
%         else
%             NLYCbCr(i,j,2) = YCbCr(i,j,2);
%             NLYCbCr(i,j,3) = YCbCr(i,j,3);
%         end
%     end
% end
% end

