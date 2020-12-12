%% TEST FACETRIANGLE ON DB1
clear;
pathBegin = 'data/DB1/db1_';
pathEnd = '.jpg';
picIndexString = '';
for i = 1:16
    if i < 10
       picIndexString = ['0' int2str(i)];
    else
       picIndexString = int2str(i);
    end
    pathString = [pathBegin picIndexString pathEnd];
    IM = imread(pathString);
    [leftEyeCoord, rightEyeCoord, mouthCoord] = faceTriangle(IM);
    subplot(4,4,i);
    imshow(IM);
    hold on;
    plot(mouthCoord(1), mouthCoord(2), '*');
    plot(leftEyeCoord(1), leftEyeCoord(2), '*');
    plot(rightEyeCoord(1), rightEyeCoord(2), '*');
end

%% TEST FACETRIANGLE ON DB2
load db2paths.mat;

picIndex = 0;
for i = 1:size(paths, 2)
    picIndex = picIndex + 1;
    if picIndex == 16
        picIndex = 1;
        figure;
    end
    currentPath = paths(i);
    IM = imread(currentPath);
    eyeM = eyeMap(IM);
    [leftEyeCoord, rightEyeCoord, mouthCoord] = faceTriangle(IM);
    subplot(4,4,picIndex);
    imshow(eyeM);
%     hold on;
%     plot(mouthCoord(1), mouthCoord(2), '*');
%     plot(leftEyeCoord(1), leftEyeCoord(2), '*');
%     plot(rightEyeCoord(1), rightEyeCoord(2), '*');
end

%% TEST FACE NORMALIZATION ON DB2
load db2paths.mat;

picIndex = 0;
for i = 1:size(paths, 2)
    picIndex = picIndex + 1;
    if picIndex == 16
        picIndex = 1;
        figure;
    end
    currentPath = paths(i);
    IM = imread(currentPath);
    normFace = normalizeFace(IM);
    faceMask = getFaceMask(IM);
    eyeM = eyeMap(IM);
    [leftEyeCoord, rightEyeCoord, mouthCoord] = faceTriangle(IM);
    subplot(4,4,picIndex);
    imshow(IM);
    hold on;
    plot(mouthCoord(1), mouthCoord(2), '*');
    plot(leftEyeCoord(1), leftEyeCoord(2), '*');
    plot(rightEyeCoord(1), rightEyeCoord(2), '*');
end

%% TEST FACE MASK ON DB2
load db2paths.mat;

picIndex = 0;
for i = 1:size(paths, 2)
    picIndex = picIndex + 1;
    if picIndex == 16
        picIndex = 1;
        figure;
    end
    currentPath = paths(i);
    IM = imread(currentPath);
    faceMask = getFaceMask(IM);
    
    subplot(4,4,picIndex);
    imshow(im2double(IM) .* faceMask);
end

%% TEST NLYCBCR
clear;
% IM = imread('data/DB2/bl_02.jpg');

load db2paths.mat;

picIndex = 0;
allMouths = [0 0];
for i = 1:size(paths, 2)
    picIndex = picIndex + 1;
    if picIndex == 16
        picIndex = 1;
        figure;
    end
    currentPath = paths(i);
    IM = imread(currentPath);
    gwIM = grayworldcorrection(IM);
    
    YCbCr = rgb2ycbcr(gwIM);
    nl_YCbCr = nlYCbCr(YCbCr);
    
    subplot(4,4,picIndex);
    imshow(horzcat(YCbCr, nl_YCbCr));
end
%% Testing pupil finding
clear;
load db2paths.mat;
picIndex = 0;
for i = 1:size(paths, 2)
    picIndex = picIndex + 1;
    if picIndex == 16
        picIndex = 1;
        figure;
    end
    currentPath = paths(i);
    IM = imread(currentPath);
    [leftEyeCoord, rightEyeCoord, mouthCoord] = faceTriangle(IM);
    gwIM = colorCorrection(IM);
    [imHeight, imWidth, ~] = size(IM);

    % Check if faceTriangle already found pupils
    R = 3;
    leftEyeMask = zeros(imHeight, imWidth);
    rightEyeMask = zeros(imHeight, imWidth);

    leftEyeMask(round(leftEyeCoord(2)), round(leftEyeCoord(1))) = 1;
    rightEyeMask(round(rightEyeCoord(2)), round(rightEyeCoord(1))) = 1;
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
        disp('Pupil found on first try');
    else
        disp('Searching for left eye pupil...');
    end
    if sum(rightEyeGrayTop(:) == 1) > 0
        rightPupilFound = 1;
        disp('Pupil found on first try');
    else
        disp('Searching for right eye pupil...');
    end

    if leftPupilFound && rightPupilFound
        leftEye = leftEyeCoord;
        rightEye = rightEyeCoord;
%         return;
    end

    hsv = rgb2hsv(im2double(gwIM));
    hue = hsv(:,:,1);
    hueBlue = hue > 200/360 & hue < 270/360; % Blue
    hueBrown = hue > 5/360 & hue < 55/360; % Brown

    hueCombined = (hueBlue + hueBrown) >= 1;

    leftEyeMask = zeros(imHeight, imWidth);
    rightEyeMask = zeros(imHeight, imWidth);

    R = 13;
    leftEyeMask(round(leftEyeCoord(2)), round(leftEyeCoord(1))) = 1;
    rightEyeMask(round(rightEyeCoord(2)), round(rightEyeCoord(1))) = 1;
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
    subplot(4,4,picIndex);
    imshow(IM)
    hold on;
    if leftPupilFound == 0
        plot(leftEye(1), leftEye(2), '*');
    end
    if rightPupilFound == 0
    plot(rightEye(1), rightEye(2), '*');
    end
    plot(leftEyeCoord(1), leftEyeCoord(2), 'o');
    plot(rightEyeCoord(1), rightEyeCoord(2), 'o');

end

%%
clear;
IM = imread('data/DB2/bl_05.jpg');
[leftEyeCoord, rightEyeCoord, mouthCoord] = faceTriangle(IM);
gwIM = colorCorrection(IM);
[imHeight, imWidth, ~] = size(IM);

% Check if faceTriangle already found pupils
R = 3;
leftEyeMask = zeros(imHeight, imWidth);
rightEyeMask = zeros(imHeight, imWidth);

leftEyeMask(round(leftEyeCoord(2)), round(leftEyeCoord(1))) = 1;
rightEyeMask(round(rightEyeCoord(2)), round(rightEyeCoord(1))) = 1;
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
    disp('Pupil found on first try');
else
    disp('Searching for left eye pupil...');
end
if sum(rightEyeGrayTop(:) == 1) > 0
    rightPupilFound = 1;
    disp('Pupil found on first try');
else
    disp('Searching for right eye pupil...');
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
leftEyeMask(round(leftEyeCoord(2)), round(leftEyeCoord(1))) = 1;
rightEyeMask(round(rightEyeCoord(2)), round(rightEyeCoord(1))) = 1;
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

imshow(IM)
hold on;
if leftPupilFound == 0
    plot(leftEye(1), leftEye(2), '*');
end
if rightPupilFound == 0
plot(rightEye(1), rightEye(2), '*');
end
plot(leftEyeCoord(1), leftEyeCoord(2), 'o');
plot(rightEyeCoord(1), rightEyeCoord(2), 'o');

%%
IM = imread('data/DB2/bl_05.jpg');
gwIM = grayworldcorrection(IM);    
YCbCr = rgb2ycbcr(gwIM);
Y = YCbCr(:,:,1);
Cb = YCbCr(:,:,2);
Cr = YCbCr(:,:,3);

% Constants attained from Face Detection in Color Images
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

% Matrixes representing low/high
Yl = uint8(Y < Kl);
Yh = uint8(Y > Kh);
Yi = uint8(Y >= Kl) .* uint8(Y <= Kh);
% calc center values Cb, Equation 7
Ba = 108;
Bb = 118;
CbCenterl = Ba + (Kl - Y) * (Bb - Ba) / (Kl - Ymin);
CbCenterh = Ba + (Y - Kh) * (Bb - Ba) / (Ymax - Kh);

% % % [CbCenterl, CbCenterh] = Ci_centers(Y);
% % % [Cb_center, Cr_center] = Ci_centers(Y);

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

% % % % Calculate cluster spread
% % % [SpreadBl, SpreadBh] = Ci_clusters(Y);
% % % SpreadB = SpreadBl .* Yl + SpreadBh .* Yh;
% % % 
% % % % calculate spread of r equation 6
% % % SpreadRl = clusterSpreadL(WLcr, Y, Ymin, Wcr, Kl);
% % % SpreadRh = clusterSpreadH(WHcr, Y, Ymax, Wcr, Kh); 
% % % SpreadR = SpreadRl .* Yl + SpreadRh .* Yh;