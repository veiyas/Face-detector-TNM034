function [combEyeMap] = eyeMap(img)

%convert to ycbcr color space
imgYcbcr = rgb2ycbcr(img);
Y = im2double(imgYcbcr(:,:,1));
Cb = im2double(imgYcbcr(:,:,2));
Cr = im2double(imgYcbcr(:,:,3));

%Constructing EyeMapC (Chrominance components)
%Cb^2
Cb2 = Cb.^2;

%Cr^2 (C with ~ sign above) which is the negative of Cr (i.e 255-Cr)
nCr = 255 - Cr;
nCr2 = nCr.^2;

%Cb/Cr
CbCr = Cb ./ Cr;

%Normalize variables to the range [0,255]
normCb2 = Cb2 - min(Cb2(:));
normCb2 = normCb2 ./ max(normCb2(:)); %Normalized between [0,1]
Cb2 = normCb2.*255; %Scale to range between [0,255]

normNCr2 = nCr2 - min(nCr2(:));
normNCr2 = normNCr2./ max(normNCr2(:));
nCr2 = normNCr2.*255; 

normCbCr = CbCr - min(CbCr(:));
normCbCr = normCbCr./ max(normCbCr(:));
CbCr = normCbCr.*255; 

eyeMapC = (1/3) * (Cb2 + nCr2 + CbCr);

%Enhance chroma eye map by histogram equalization
eyeMapC = histeq(uint8(eyeMapC));

%Constructing EyeMapL (Luminance components)
%We can use morphological operators such as dilation and erosion to
%emphasize brighter and darker pixels in the luma component around eyes.

%Construct a structuring element for dilation and erosion
radius1 = 4;
SE1 = strel('disk', radius1);

dilation = imdilate(Y,SE1);
erosion = imerode(Y,SE1);

eyeMapL = dilation./(erosion+1);

%Combines eyeMap
combEyeMap = imfuse(eyeMapL, eyeMapC, 'blend');

SE2 = strel('disk', 5);
combEyeMap = imdilate(combEyeMap, SE2);



