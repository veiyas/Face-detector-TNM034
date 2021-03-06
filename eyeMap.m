function [eyeMap] = eyeMap(img)
% Return binary image
YCbCr = rgb2ycbcr(img);
faceMask = getFaceMask(img);

Y = im2double(YCbCr(:,:,1));
Cb = im2double(YCbCr(:,:,2));
Cr = im2double(YCbCr(:,:,3));

Cbsqr = Cb.^2;
CbNegsqr = 1 - Cr.^2;
CbCr = Cb ./ Cr;

eyeMapC = (1/3)*(Cbsqr + CbNegsqr + CbCr);
eyeMapC = histeq(eyeMapC);

diskSize = 7;
kernel = strel('disk', diskSize);

nominator = imdilate(Y, kernel);
denominator = imerode(Y, kernel);

eyeMapL = nominator ./ (denominator + 1);

eyeMap = eyeMapC .* eyeMapL;
eyeMap = imdilate(eyeMap, kernel);
eyeMap = eyeMap .* faceMask;
eyeMap = im2uint8(eyeMap);
eyeMapStretch = imadjust(eyeMap, stretchlim(eyeMap), []);
eyeMapBin = eyeMapStretch > 250;

closer = strel('disk', 5);
denoiser = strel('disk', 4);
vertLinesRemover = strel('line', 10, 0);

eyeMapBinMorphed = imopen(eyeMapBin, denoiser);
eyeMapBinMorphed = imopen(eyeMapBinMorphed, vertLinesRemover);
eyeMapBinMorphed = imclose(eyeMapBinMorphed, closer);
eyeMap = eyeMapBinMorphed;

% %convert to ycbcr color space
% imgYcbcr = rgb2ycbcr(img);
% Y = im2double(imgYcbcr(:,:,1));
% Cb = im2double(imgYcbcr(:,:,2));
% Cr = im2double(imgYcbcr(:,:,3));
% 
% %Constructing EyeMapC (Chrominance components)
% %Cb^2
% Cb2 = Cb.^2;
% 
% %Cr^2 (C with ~ sign above) which is the negative of Cr (i.e 255-Cr)
% nCr = 255 - Cr;
% nCr2 = nCr.^2;
% 
% %Cb/Cr
% CbCr = Cb ./ Cr;
% 
% %Normalize variables to the range [0,255]
% normCb2 = Cb2 - min(Cb2(:));
% normCb2 = normCb2 ./ max(normCb2(:)); %Normalized between [0,1]
% Cb2 = normCb2.*255; %Scale to range between [0,255]
% 
% normNCr2 = nCr2 - min(nCr2(:));
% normNCr2 = normNCr2./ max(normNCr2(:));
% nCr2 = normNCr2.*255; 
% 
% normCbCr = CbCr - min(CbCr(:));
% normCbCr = normCbCr./ max(normCbCr(:));
% CbCr = normCbCr.*255; 
% 
% eyeMapC = (1/3) * (Cb2 + nCr2 + CbCr);
% 
% %Enhance chroma eye map by histogram equalization
% eyeMapC = histeq(eyeMapC);
% 
% %Constructing EyeMapL (Luminance components)
% %We can use morphological operators such as dilation and erosion to
% %emphasize brighter and darker pixels in the luma component around eyes.
% 
% %Construct a structuring element for dilation and erosion
% size = 8;
% SE1 = strel('disk', size);
% 
% dilation = imdilate(Y,SE1);
% erosion = imerode(Y,SE1);
% 
% eyeMapL = dilation./(erosion+1);
% 
% %Combines eyeMap
% %combEyeMap = imfuse(eyeMapL, eyeMapC, 'blend');
% combEyeMap = eyeMapL.*eyeMapC;
% 
% % SE2 = strel('disk', 3);
% % combEyeMap = imdilate(combEyeMap, SE2);
% 
% combEyeMap = uint8(rescale(combEyeMap,0,255));
end



