function [leftEye, rightEye] = findPupils(IM, leftEyeCoord, rightEyeCoord)
%findPupils Do an extra search for the pupils around the parameter points
% This function checks for pixel values typically found in a small area around the eye
gwIM = colorCorrection(IM);
[imHeight, imWidth, ~] = size(IM);

% Check if faceTriangle already found pupils
R = 3;
leftEyeMask = zeros(imHeight, imWidth);
rightEyeMask = zeros(imHeight, imWidth);

leftEyeMask = bwdist(leftEyeMask) <= R;
rightEyeMask = bwdist(rightEyeMask) <= R;

leftEyeMaskedIm = im2double(im2double(gwIM) .* leftEyeMask);
rightEyeMaskedIm = im2double(im2double(gwIM) .* rightEyeMask);

leftEyeGrayTop = rgb2gray(leftEyeMaskedIm) > 0.9;
rightEyeGrayTop = rgb2gray(rightEyeMaskedIm) > 0.9;

leftPupilFound = 0;
rightPupilFound = 0;

if sum(leftEyeGrayTop(:) == 1) > 0
    leftPupilFound = 1;
end
if sum(rightEyeGrayTop(:) == 1) > 0
    rightPupilFound = 1;
end

if leftPupilFound && rightPupilFound
    leftEye = leftEyeCoord;
    rightEye = rightEyeCoord;
    return;
end

hsv = rgb2hsv(im2double(gwIM));
hue = hsv(:,:,1);
hueBlue = hue > 200/360 & hue < 270/360; % Blue
hueBrown = hue > 5/360 & hue < 55/360; % Brown

hueCombined = (hueBlue + hueBrown) >= 1;

leftEyeMask = zeros(imHeight, imWidth);
rightEyeMask = zeros(imHeight, imWidth);

R = 13;
leftEyeMask = bwdist(leftEyeMask, 'quasi-euclidean') <= R;
rightEyeMask = bwdist(rightEyeMask, 'quasi-euclidean') <= R;


leftEyeMaskedIm = im2double(hueCombined .* leftEyeMask);
rightEyeMaskedIm = im2double(hueCombined .* rightEyeMask);

onlyLeftEye = imfill(leftEyeMaskedIm,'holes') .* im2double(gwIM);
onlyRightEye = imfill(rightEyeMaskedIm,'holes') .* im2double(gwIM);

leftEyeGray = rgb2gray(onlyLeftEye);
leftEyeGray = contrastStretch(leftEyeGray);
rightEyeGray = rgb2gray(onlyRightEye);
rightEyeGray = contrastStretch(rightEyeGray);

% subplot(1,2,1);
% imshow(leftEyeGray > 0.999);
% hold on
% plot(leftEyeCoord(1), leftEyeCoord(2), 'o');
% plot(rightEyeCoord(1), rightEyeCoord(2), 'o');
% subplot(1,2,2);
% imshow(rightEyeGray > 0.999);
% hold on
% plot(leftEyeCoord(1), leftEyeCoord(2), 'o');
% plot(rightEyeCoord(1), rightEyeCoord(2), 'o');

leftPupil = leftEyeGray > 0.999;
rightPupil = rightEyeGray > 0.999;

% YCbCrLEFT = rgb2ycbcr(onlyLeftEye);
% YCbCrRIGHT = rgb2ycbcr(onlyRightEye);
% CrLEFT = YCbCrLEFT(:,:,3).^2;
% CrRIGHT = YCbCrRIGHT(:,:,3).^2;
% 
% CrStrLeft = imadjust(CrLEFT, stretchlim(CrLEFT));
% CrStrRight = imadjust(CrRIGHT, stretchlim(CrRIGHT));
% leftPupil = CrStrLeft < 0.24;
% rightPupil = CrStrRight < 0.24;

leftEyeStats = regionprops('table', leftPupil, 'Centroid');
leftEyeCentroids = cat(1,leftEyeStats.Centroid);

rightEyeStats = regionprops('table', rightPupil, 'Centroid');
rightEyeCentroids = cat(1,rightEyeStats.Centroid);

leftEyeIsSet = leftPupilFound;
rightEyeIsSet = rightPupilFound;
% If we dont find any good eye regions we use faceTriangle values
if (size(leftEyeCentroids, 1) == 0 && leftEyeIsSet == 0) || (size(rightEyeCentroids, 1) == 0 && rightEyeIsSet == 0)
    if size(leftEyeCentroids, 1) == 0
        leftEye = leftEyeCoord;
        leftEyeIsSet = 1;
    end
    if size(rightEyeCentroids, 1) == 0
        rightEye = rightEyeCoord;
        rightEyeIsSet = 1;
    end
end
if (size(leftEyeCentroids, 1) == 1  && leftEyeIsSet == 0) && (size(rightEyeCentroids, 1) == 1 && rightEyeIsSet == 0)
    leftEye = leftEyeCentroids(1,:);
    rightEye = rightEyeCentroids(1,:);
else
    if leftEyeIsSet == 0
        avgCentroid = [0 0];
        for i = 1:size(leftEyeCentroids,1)
            avgCentroid = avgCentroid + leftEyeCentroids(i,:);
        end
        leftEye = avgCentroid ./ i;
    end
    
    if rightEyeIsSet == 0
        avgCentroid = [0 0];
        for i = 1:size(rightEyeCentroids,1)
            avgCentroid = avgCentroid + rightEyeCentroids(i,:);
        end
        rightEye = avgCentroid ./ i;
    end
end

% subplot(4,4,picIndex);
% imshow(IM)
% hold on;
if leftPupilFound == 1
    leftEye = leftEyeCoord;
end
if rightPupilFound == 1
    rightEye = rightEyeCoord;
end
% plot(leftEyeCoord(1), leftEyeCoord(2), 'o');
% plot(rightEyeCoord(1), rightEyeCoord(2), 'o');
end

