function [faceMask] = getFaceMask(IM)
% Calculate lighting corrected image (gray world assumption)
alpha = mean(IM(:,:,2), 'all') / mean(IM(:,:,1), 'all');
beta = mean(IM(:,:,2), 'all') / mean(IM(:,:,3), 'all');

gwIM(:,:,1) = alpha .* IM(:,:,1);
gwIM(:,:,2) = IM(:,:,2);
gwIM(:,:,3) = beta .* IM(:,:,3);
imSize = size(gwIM);

% Determine if lighting correction is needed
grayIM = rgb2gray(IM);
% maxVal = max(grayIM(:));
% top95binIM = grayIM >= maxVal*0.95;
% nTop95vals = sum(top95binIM(:) == 1);
% totPixels = imSize(1)*imSize(2);


YCbCr = rgb2ycbcr(gwIM);
nl_YCbCr = nlYCbCr(YCbCr);

theta = 2.53;
cx = 109.38;
cy = 152.02;
ecx = 1.6;
ecy = 2.41;

imgCb(:,:) = double(nl_YCbCr(:,:,2)) - cx;
imgCr(:,:) = double(nl_YCbCr(:,:,3)) - cy;

posX = imgCb*cos(theta)+imgCr*sin(theta);
posY = imgCb*-sin(theta)+imgCr*cos(theta);

% Check if points are inside ellipse
a2 = 25.39 * 25.39;
b2 = 14.03*14.03;

threshold = 2.0;
faceMask = ( ((posX - ecx).^2 ./ a2 + (posY - ecy).^2 ./ b2) <= threshold);

% morphological operations, closing, open, closing
diskSize = 20;
kernel = strel('disk', diskSize);
faceMask = imclose(faceMask, kernel);

diskSize = 4;
kernel = strel('disk', diskSize);
faceMask = imopen(faceMask, kernel);

diskSize = 25;
kernel = strel('disk', diskSize);
faceMask = imclose(faceMask, kernel);

% fill in blanks of the objects
faceMask = imfill(faceMask,'holes');

largest = bwareafilt(faceMask,1);

% assume largest area is face and use only that using matlab magic
% if (sum(sum(faceMask)) == 0)
%     faceMask = skinModel(transformedImage, threshold+.2);
%     return
% else
%     faceMask = largestArea(faceMask);
%     faceMask(:,:,2) = faceMask;
%     faceMask(:,:,3) = faceMask(:,:,2);
% end
end

