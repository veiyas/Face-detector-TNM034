function bbRect = findFaceBoundingBox(rgbImage)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
IM = rgbImage;

% Calculate lighting corrected image (gray world assumption)
alpha = mean(IM(:,:,2), 'all') / mean(IM(:,:,1), 'all');
beta = mean(IM(:,:,2), 'all') / mean(IM(:,:,3), 'all');

gwIM(:,:,1) = alpha .* IM(:,:,1);
gwIM(:,:,2) = IM(:,:,2);
gwIM(:,:,3) = beta .* IM(:,:,3);
imSize = size(gwIM);

% Determine if lighting correction is needed
grayIM = rgb2gray(IM);
maxVal = max(grayIM(:));
top95binIM = grayIM >= maxVal*0.95;
nTop95vals = sum(top95binIM(:) == 1);

totPixels = imSize(1)*imSize(2);
correctionTolerance = 0.3; % How much of image needs to be bright to corrected
if (nTop95vals/totPixels) < correctionTolerance
    gwIM = IM;
end

% Construct skin color boundary in YCbCr space
% Not used atm, can perhaps use the elliptical model instead
skinToneCB = [85 100 110 130 122 108 108 98 85];
skinToneCR = [160 140 125 138 150 158 162 175 160];

YCbCr=rgb2ycbcr(gwIM);
% DO NON LINEAR TRANSFORMATION OF THE YCBCR COLOR SPACE HERE

% Construct mask
mask = zeros(imSize(1), imSize(2));
for j = 1:imSize(1)
    for k = 1:imSize(2)
        tempCb = double(YCbCr(j,k,2));
        tempCr = double(YCbCr(j,k,3));        
        
%         if inpolygon(tempCb, tempCr, skinToneCB, skinToneCR)
%             mask(j,k) = 1;
%         end
        if (tempCr > 130 && tempCr < 160 && tempCb > 80 && tempCb < 145)
            mask(j,k) = 1;
        end
    end
end

% Attempt to remove noise and background stuff
% Can probably be expanded
dotRemover = strel('disk', 4);
tallObjectRemover = strel('line', 50, 0);
wideObjectRemover = strel('line', 50, 90);
openedMask = imopen(mask, dotRemover);
openedMask = imopen(openedMask, tallObjectRemover);
openedMask = imopen(openedMask, wideObjectRemover);

% Get image stats
stats = regionprops(openedMask, 'BoundingBox', 'EulerNumber');
bbRect = stats.BoundingBox;

% Check if bounding box is too small (failure), set a big bounding box
% This is a hack to aid testing
if (bbRect(3)*bbRect(4)) / totPixels < 0.3 || (bbRect(3)*bbRect(4)) / totPixels > 0.9
   centerX = round(imSize(2) / 2);
   centerY = round(imSize(1) / 2);
   lengthX = centerX / 2;
   lengthY = centerY / 2;
   bbRect = [centerX - lengthX centerY - lengthY centerX centerY];
end

faceMask = bbRect;

end

