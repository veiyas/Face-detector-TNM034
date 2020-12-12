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
    grayIM = im2double(rgb2gray(IM));
    faceMask = getFaceMask(IM);
    
    subplot(4,4,picIndex);
    imshow(grayIM .* faceMask);
end

%% Testing pupil finding
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
    [leftEyeCoord, rightEyeCoord, mouthCoord] = faceTriangle(IM);
%     subplot(4,4,picIndex);
%     imshow(IM)
%     hold on;
% %     plot(leftEye(1), leftEye(2), '*');
% %     plot(rightEye(1), rightEye(2), '*');
%     plot(leftEyeCoord(1), leftEyeCoord(2), 'o');
%     plot(rightEyeCoord(1), rightEyeCoord(2), 'o');
%     plot(mouthCoord(1), mouthCoord(2), 'o');
    allMouths = allMouths + mouthCoord;
end
avgMouth = allMouth ./ i
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
% Get image size
[imgHeight, imgWidth] = size(gwIM);

% Get mouth map and mouth coords
mouthmap = mouthMap(gwIM);
mouthCoord = getMouthCoord(mouthmap);

eyeMapBinMorphed = eyeMap(gwIM);

eyeMapMouthFix = eyeMapBinMorphed;
% 
% if isempty(mouthCoord)
% else

eyeMapMouthFix(round(mouthCoord(2))-30:imgHeight, :) = 0;


% Get all potential eye candidates
eyeStats = regionprops('table', eyeMapMouthFix, 'Centroid', 'BoundingBox', 'Area');
% Store all centroids in a vector
eyeCentroids = cat(1,eyeStats.Centroid);
eyeAreas = cat(1, eyeStats.Area);
eyesDetected = size(eyeCentroids);

% Remove all eyes that does not fulfill certain attributes
ple = []; % Potential left eyes
pleArea = []; % Potenatial left eye areas
pleAngle = [];
l = 1; % Left eye counter
pre = []; % Potential right eyes
preArea = []; % Potenatial right eye areas
preAngle = [];
r = 1; % Right eye counter

heightDiffTolerance = 30;
areaDiffTolerance = 670;
normDiffTolerance = 20;
angleDiffTolerance = 18;
smallestEyeMouthAngle = 12;
largestEyeMouthAngle = 45;
mouthNormal = [0 1];

% Filter out candidates that are:
%   Not on correct side of mouth
%   Too far apart vertically
%   Not in similar size
%   Vector to mouth should be roughly same length
%   Vector to mouth has similar angles to [0 1]
%   Vector to mouth has reasonably large angle
for i = 1:eyesDetected-1
    for j = i+1:eyesDetected
        tmpEyeArea1 = eyeAreas(i);
        tmpEyeArea2 = eyeAreas(j);
        
        tmpEyeMtrx = [eyeCentroids(i,:); eyeCentroids(j,:)];
        
        [~,leftEyeIndex] = min(tmpEyeMtrx, [], 1);
        [~,rightEyeIndex] = max(tmpEyeMtrx, [], 1);
        
        tmpLeftEye = tmpEyeMtrx(leftEyeIndex(1),:);
        tmpRightEye = tmpEyeMtrx(rightEyeIndex(1),:);
        
        leftEyeMouthVec = mouthCoord - tmpLeftEye;
        leftEyeMouthNorm = norm(leftEyeMouthVec);
        leftEyeMouthAngle = acosd(dot(leftEyeMouthVec, mouthNormal) /leftEyeMouthNorm);
        
        rightEyeMouthVec = mouthCoord - tmpRightEye;
        rightEyeMouthNorm = norm(rightEyeMouthVec);
        rightEyeMouthAngle = acosd(dot(rightEyeMouthVec, mouthNormal) /rightEyeMouthNorm);
        
        % Filters
        correctSides = tmpLeftEye(1) < mouthCoord(1) && mouthCoord(1) < tmpRightEye(1);
        belowMouth = tmpLeftEye(2) < mouthCoord(2) && tmpRightEye(2) < mouthCoord(2);
        smallVertDiff = abs(tmpLeftEye(1,2) - tmpRightEye(1,2)) < heightDiffTolerance;
        smallAreaDiff = abs(tmpEyeArea1 - tmpEyeArea2) < areaDiffTolerance;
        smallMouthNormDiff = abs(leftEyeMouthNorm - rightEyeMouthNorm) < normDiffTolerance;
        smallEyeMouthAngleDiff = abs(leftEyeMouthAngle - rightEyeMouthAngle) < angleDiffTolerance;
        reasonableAngle = leftEyeMouthAngle > smallestEyeMouthAngle && rightEyeMouthAngle > smallestEyeMouthAngle;
        reasonableAngle = reasonableAngle && leftEyeMouthAngle < largestEyeMouthAngle && rightEyeMouthAngle < largestEyeMouthAngle;
        
        if correctSides && belowMouth && smallVertDiff && smallMouthNormDiff && smallEyeMouthAngleDiff && smallAreaDiff && reasonableAngle
%             leftEyeMouthVec
%             rightEyeMouthVec
%             leftEyeMouthAngle
%             rightEyeMouthAngle
            ple(l,:) = tmpLeftEye;
            pleAreas(l) = tmpEyeArea1;
            pleAngle(l) = leftEyeMouthAngle;
            l = l+1;
            pre(r,:) = tmpRightEye;
            preAreas(r) = tmpEyeArea2;
            preAngle(r) = rightEyeMouthAngle;
            r = r+1;
        end
    end
end
% If we only have 1 candidate left its probably correct
if size(ple,1) == 1 && size(pre,1) == 1
    leftEyeCoord = ple;
    rightEyeCoord = pre;
elseif size(ple,1) == 0 && size(pre,1) == 0 % Guesstimated triangle, is this allowed?
    disp('Could not find eyes');
    leftEyeCoord = mouthCoord - [60 120];
    rightEyeCoord = mouthCoord - [-60 120];
else % Check pair by pair for the lowest angle to mouthNormal
    lowestAngle = 10000;
    lowestAngleIndex = -1;
    
    for i = 1:size(pleAngle,2)
        tmpLeftEyeAngle = pleAngle(i);
        tmpRightEyeAngle = preAngle(i);
        
        if min(tmpLeftEyeAngle,tmpRightEyeAngle) < lowestAngle
            lowestAngle = min(tmpLeftEyeAngle,tmpRightEyeAngle);
            lowestAngleIndex = i;
        end
    end
    pleAngle(lowestAngleIndex)
    preAngle(lowestAngleIndex)
    
    leftEyeCoord = ple(lowestAngleIndex,:);
    rightEyeCoord = pre(lowestAngleIndex,:);
end


subplot(1,2,1);
imshow(eyeMapMouthFix);
hold on;
plot(mouthCoord(1), mouthCoord(2), 'x');
plot(mouthCoord(1), mouthCoord(2), '*');
plot(leftEyeCoord(1), leftEyeCoord(2), '*');
plot(rightEyeCoord(1), rightEyeCoord(2), '*');
subplot(1,2,2);
imshow(IM);
hold on;
plot(mouthCoord(1), mouthCoord(2), '*');
plot(leftEyeCoord(1), leftEyeCoord(2), '*');
plot(rightEyeCoord(1), rightEyeCoord(2), '*');



%%

for j=1:eyesDetected
    eye = eyeCentroids(j,:);
    % If left side of mouth and within resonable distance.
    if(eye(1) < mouthCoord(1) && abs(eye(1)-mouthCoord(1)) < 75)
        ple(l,:) = eye;
        l = l+1;
    % If instead right side of mouth and within resonable distance.  
    elseif (eye(1) > mouthCoord(1) && abs(eye(1)-mouthCoord(1)) < 75)
        pre(r,:) = eye;
        r = r+1;
    end
end

% Choose the eye pair with the least y-coord difference
bestPair = [];
yDiff = 1000;
for j=1:l
    for k=1:r
        if(size(ple,1) > 0 && size(pre,1) > 0)
            if(abs(ple(j,2) - pre(k,2)) < yDiff)
                bestPair = cat(2,ple(j),pre(k));
            end
        end
    end
end
combo(i,:) = cat(2,mouth,bestPair);    





% Somehow choose the best potential combo and return eyes and mouth

%temporary solution. Just pick the first one.
if(size(combo,1) > 0)
    mouthPos = combo(1,1:2);
    leftEye = combo(1,3:4);
    rightEye = combo(1,5:6);
else
    leftEye = [0 0];
    rightEye = [0 0];
    mouthPos = [0 0];
end



%%
% Attempt to find eyes
s = regionprops(eyeMapMouthFix, 'centroid', 'Circularity');
eyeCentroids = cat(1, s.Centroid);
eyeCircs = cat(1, s.Circularity);

nRegions = size(eyeCircs);

eyeHeightDiffTolerance = 20;
% Eyes should be at roughly same height and round-ish
for i = 1:nRegions-1
    for j = i+1 : nRegions
        if eyeCentroids(i,2) - eyeCentroids(j,2) < eyeHeightDiffTolerance
            
        end
    end
end


subplot(2,2,1);
imshow(IM);
hold on
plot(mouthCoord(1), mouthCoord(2), 'x');
hold off
subplot(2,2,2);
imshow(eyeMapStretch);
subplot(2,2,3);
imshow(eyeMapMouthFix);




% imshow(im2double(IM) .* faceMask)
% imshow(eyemap)


%%
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

% totPixels = imSize(1)*imSize(2);
% correctionTolerance = 0.3; % How much of image needs to be bright to corrected
% if (nTop95vals/totPixels) < correctionTolerance
%     gwIM = IM;
% end

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

threshold = 1.8;
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

bwIM = im2double(rgb2gray(IM));
maskedIM = im2double(gwIM) .* faceMask;
eyemap = eyeMap(gwIM);
eyemap = eyemap .* faceMask;
figure
% imshow(bwIM .* largest);
imshow(eyemap);

%%
clear;
IM = imread('data/DB2/il_01.jpg');
% IM = imread('data/DB2/cl_10.jpg');

% Calculate lighting corrected image (gray world assumption)
alphaGW = mean(IM(:,:,2), 'all') / mean(IM(:,:,1), 'all');
betaGW = mean(IM(:,:,2), 'all') / mean(IM(:,:,3), 'all');

gwIM(:,:,1) = alphaGW .* IM(:,:,1);
gwIM(:,:,2) = IM(:,:,2);
gwIM(:,:,3) = betaGW .* IM(:,:,3);
imSize = size(gwIM);

% White patch
alphaWP = max(max(IM(:,:,2))) / max(max(IM(:,:,1)));
betaWP = max(max(IM(:,:,2))) / max(max(IM(:,:,3)));

wpIM(:,:,1) = alphaWP .* IM(:,:,1);
wpIM(:,:,2) = IM(:,:,2);
wpIM(:,:,3) = betaWP .* IM(:,:,3);

% subplot(1,3,1)
% imshow(IM);
% subplot(1,3,2)
% imshow(gwIM);
% subplot(1,3,3)
% imshow(wpIM);

% Determine if lighting correction is needed
grayIM = rgb2gray(IM);
maxVal = max(grayIM(:));
top95binIM = grayIM >= maxVal*0.95;
nTop95vals = sum(top95binIM(:) == 1);

totPixels = imSize(1)*imSize(2);
correctionTolerance = 0.3; % How much of image needs to be bright to corrected
if (nTop95vals/totPixels) < correctionTolerance
    gwIM = IM;
    disp('Image was lighting corrected')
end

AInv = imcomplement(gwIM);
[BInv, TInv] = imreducehaze(AInv, 'Method', 'approxdcp', 'ContrastEnhancement', 'none');
lightCorr = TInv .* im2double(grayIM);
corrIM = imadjust(lightCorr,stretchlim(lightCorr),[]);
gwIM = imcomplement(BInv);

%%

% Construct skin color boundary in YCbCr space
% Not used atm, can perhaps use the elliptical model instead
skinToneCBorig = [85 100 110 130 122 108 108 98 85];
skinToneCRorig = [160 140 125 138 150 158 162 175 160];

p = polyshape(skinToneCBorig, skinToneCRorig);
[xc_P, yc_P] = centroid(p);
shrinkFactor = 0.5;
scaledPoly = scale(p, shrinkFactor, [xc_P yc_P]);

YCbCr=rgb2ycbcr(gwIM);
Kl = 125;
Kh = 188;
WCb = 46.97;
WCr = 38.76;

NLYCbCr = zeros(imSize(1), imSize(2), imSize(3));
NLYCbCr = uint8(NLYCbCr);
NLYCbCr(:,:,1) = YCbCr(:,:,1);
for i = 1:imSize(1)
    for j = 1:imSize(2)
        Y = YCbCr(i,j,1);
        
        if Y < Kl || Kh < Y
            [CbC, CrC] = Ci_centers(Y); % Centers
            [CbCK, CrCK] = Ci_centers(Kh); % Centers with Kh as input
            [CbW, CrW] = Ci_weights(Y); % Weights
            NLCb = (Y - CbC)*(WCb / CbW) + CbCK;
            NLCr = (Y - CrC)*(WCr / CrW) + CrCK;
            NLYCbCr(i,j,2) = uint8(NLCb);
            NLYCbCr(i,j,3) = uint8(NLCr);
        else
            NLYCbCr(i,j,2) = YCbCr(i,j,2);
            NLYCbCr(i,j,3) = YCbCr(i,j,3);
        end
    end
end

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
        if (tempCr > 135 && tempCr < 150 && tempCb > 115 && tempCb < 146)
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
stats = regionprops(openedMask, 'BoundingBox', 'EulerNumber');
bbRect = stats.BoundingBox;
statsSize = size(stats, 1);
imshow(IM);
hold on;
rectangle('Position', bbRect, 'EdgeColor','r', 'LineWidth', 3)

