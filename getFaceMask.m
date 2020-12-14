function [largest] = getFaceMask(IM)
gwIM = grayworldcorrection(IM);
YCbCr = rgb2ycbcr(gwIM);
nl_YCbCr = nlYCbCr(YCbCr);

% Construct ellipse
theta = 2.53;
cx = 109.38;
cy = 152.02;
ecx = 1.6;
ecy = 2.41;

a2 = 25.39 * 25.39;
b2 = 14.03 * 14.03;

imgCb(:,:) = double(nl_YCbCr(:,:,2)) - cx;
imgCr(:,:) = double(nl_YCbCr(:,:,3)) - cy;

posX = imgCb * cos(theta) + imgCr * sin(theta);
posY = imgCb * -sin(theta) + imgCr * cos(theta);

% Check which skin samples are inside ellipse
threshold = 2.0;
faceMask = ( ((posX - ecx).^2 ./ a2 + (posY - ecy).^2 ./ b2) <= threshold);

% morphological operations
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

% We assume that the largest area left is the face
largest = bwareafilt(faceMask,1);
end

