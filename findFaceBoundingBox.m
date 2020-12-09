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
% if (nTop95vals/totPixels) < correctionTolerance
%     gwIM = IM;
% end

% Construct skin color boundary in YCbCr space
% Not used atm, can perhaps use the elliptical model instead
skinToneCBorig = [85 100 110 130 122 108 108 98 85];
skinToneCRorig = [160 140 125 138 150 158 162 175 160];

p = polyshape(skinToneCBorig, skinToneCRorig);
[xc_P, yc_P] = centroid(p);
shrinkFactor = 0.55;
scaledPoly = scale(p, shrinkFactor, [xc_P yc_P]);

YCbCr=rgb2ycbcr(gwIM);
Kl = 125;
Kh = 188;
WCb = 46.97;
WCr = 38.76;

NLYCbCr = nlYCbCr(YCbCr);

% Elliptical skin tone boundary
theta = 2.53;
t = 1:1:256;
cx = 109.38;
cy = 152.02;
cossinMatrix = [cos(theta) sin(theta); -sin(theta) cos(theta)];
CbCrVec = [t - cx; t - cy];
res = cossinMatrix * CbCrVec;

% subplot(1,2,1)
% imshow(YCbCr);
% subplot(1,2,2)
% imshow(NLYCbCr);

% Construct mask
mask = zeros(imSize(1), imSize(2));
for j = 1:imSize(1)
    for k = 1:imSize(2)
        tempCb = double(NLYCbCr(j,k,2));
        tempCr = double(NLYCbCr(j,k,3));        
        
%         if isinterior(scaledPoly, tempCb, tempCr)
%             mask(j,k) = 1;
%         end
        if (tempCr > 135 && tempCr < 150 && tempCb > 110 && tempCb < 146)
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
% openedMask = imopen(openedMask, tallObjectRemover);
% openedMask = imopen(openedMask, wideObjectRemover);

% Get image stats
stats = regionprops(openedMask, 'BoundingBox');
SizeOfFirstField = size(stats, 1);
if SizeOfFirstField == 0
    bbRect = [0 0 10 10];
else
    bbRect = stats.BoundingBox;
end

end

