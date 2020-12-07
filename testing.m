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
    grayIM = im2double(rgb2gray(IM));
    mouthmap = mouthMap(IM);
    
    %Get all potential mouth candidates
    mouthStats = regionprops('table', mouthmap, 'Centroid', 'BoundingBox', 'Area');
    mouthCentroids = cat(1, mouthStats.Centroid);
    mouthAreas = cat(1, mouthStats.Area);
    mouthsDetected = size(mouthCentroids);
    [imgHeight, imgWidth] = size(IM);

% The region with the largest area is 99% the mouth
% So we search for it with some extra filters as precaution
    mouthCoord = [];
    currArea = -1;
    k = 1;
    for j = 1:mouthsDetected
        % The mouth should be wider than its height
        if(mouthStats.BoundingBox(j,3) > mouthStats.BoundingBox(j,4))
            % The mouth should not be in the upper part of the image
            % Maybe even bottom half?
            if(mouthCentroids(j,2) > imgHeight/2)
                % The region needs to have a larger area
                if (mouthAreas(j) > currArea)
                    currArea = mouthAreas(j);
                    mouthCoord = mouthCentroids(j,:);
                    k = k+1;
                end            
            end
        end
    end
    
%     res = horzcat(grayIM, mm);
    subplot(4,4,i);
    imshow(IM);
    hold on;
    plot(mouthCoord(1), mouthCoord(2), 'x');
end

%%
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
    gwIM = grayworldcorrection(IM);
    grayIM = im2double(rgb2gray(IM));
    faceMask = getFaceMask(IM);
    
    maskedIM = im2double(IM) .* faceMask;
    
    subplot(4,4,i);
    imshow(grayIM .* faceMask);
    hold on;
    
    [left_eye, right_eye, mouth, numberOfEyes] = get_eye_mouth_coord(gwIM, faceMask);
    if numberOfEyes == 2
        plot(left_eye(1), left_eye(2), 'o');
        plot(right_eye(1), right_eye(2), 'o');
        plot(mouth(1), mouth(2), 'o');
    end    
end

%%
load db2paths.mat;

picIndex = 0;
for i = 2:size(paths, 2)
    picIndex = picIndex + 1;
    if picIndex == 16
        picIndex = 1;
        figure;
    end
    currentPath = paths(i);
    IM = imread(currentPath);
    grayIM = im2double(rgb2gray(IM));
    mouthmap = mouthMap(IM);
    
    %Get all potential mouth candidates
    mouthStats = regionprops('table', mouthmap, 'Centroid', 'BoundingBox', 'Area');
    mouthCentroids = cat(1, mouthStats.Centroid);
    mouthAreas = cat(1, mouthStats.Area);
    mouthsDetected = size(mouthCentroids);
    [imgHeight, imgWidth] = size(IM);

% The region with the largest area is 99% the mouth
% So we search for it with some extra filters as precaution
    mouthCoord = [];
    currArea = -1;
    k = 1;
    for j = 1:mouthsDetected
        % The mouth should be wider than its height
        if(mouthStats.BoundingBox(j,3) > mouthStats.BoundingBox(j,4))
            % The mouth should not be in the upper part of the image
            % Maybe even bottom half?
            if(mouthCentroids(j,2) > imgHeight/2)
                % The region needs to have a larger area
                if (mouthAreas(j) > currArea)
                    currArea = mouthAreas(j);
                    mouthCoord = mouthCentroids(j,:);
                    k = k+1;
                end            
            end
        end
    end
    
%     res = horzcat(grayIM, mm);
    subplot(4,4,picIndex);
    imshow(IM);
    hold on;
    plot(mouthCoord(1), mouthCoord(2), 'x');
end


%%
load db2paths.mat;

picIndex = 0;
for i = 2:size(paths, 2)
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


%%
clear;
IM = imread('data/DB2/bl_01.jpg');
% IM = imread('data/DB1/db1_01.jpg');
gwIM = grayworldcorrection(IM);
YCbCr = rgb2ycbcr(gwIM);
faceMask = getFaceMask(gwIM);

% Get image size
[imgHeight, imgWidth] = size(gwIM);

% Get mouth map and mouth coords
mouthmap = mouthMap(IM);
mouthCoord = getMouthCoord(mouthmap);

% Eye map stuff
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

% Remove all data below the mouth
eyeMapMouthFix = eyeMapBinMorphed;
eyeMapMouthFix(mouthCoord(2)-30:imgHeight, :) = 0;

% Get all potential eye candidates
eyeStats = regionprops('table', eyeMapMouthFix, 'Centroid', 'BoundingBox', 'Area')
% Store all centroids in a vector
eyeCentroids = cat(1,eyeStats.Centroid);
eyeAreas = cat(1, eyeStats.Area);
eyesDetected = size(eyeCentroids);

% Remove all eyes that does not fulfill certain attributes
ple = []; % Potential left eyes
l = 1; % Left eye counter
pre = []; % Potential right eyes
r = 1; % Right eye counter

heightDiffTolerance = 50;
areaDiffTolerance = 150;
combo = []; % Stores each possible eye+mouth combo

% Filter out candidates that are:
%   Not on correct side of mouth
%   Too far apart vertically
%   Not in similar size
for i = 1:eyesDetected-1
    for j = i+1:eyesDetected
%         tmpEye1 = ;
        tmpEyeArea1 = eyeAreas(i);
        
%         tmpEye2 = ;
        tmpEyeArea2 = eyeAreas(j);
        
        tmpEyeMtrx = [eyeCentroids(i,:); eyeCentroids(j,:)];
        
        [~,leftEyeIndex] = min(tmpEyeMtrx, [], 1);
        [~,rightEyeIndex] = max(tmpEyeMtrx, [], 1);
        
        tmpLeftEye = tmpEyeMtrx(leftEyeIndex(1),:);
        tmpRightEye = tmpEyeMtrx(rightEyeIndex(1),:);
        
        if tmpLeftEye(1) < mouthCoord(1) && mouthCoord(1) < tmpRightEye(1) % Right sides of mouth
            if tmpLeftEye(2) < mouthCoord(2) && tmpRightEye(2) < mouthCoord(2) % Above mouth
                if abs(tmpLeftEye(1,2) - tmpRightEye(1,2)) < heightDiffTolerance
                    if abs(tmpEyeArea1 - tmpEyeArea2) < areaDiffTolerance
                        ple(l,:) = tmpLeftEye;
                        l = l+1;
                        pre(r,:) = tmpRightEye;
                        r = r+1;
                    end
                end
            end
        end
    end
end
ple
pre

subplot(1,2,1);
imshow(IM);
hold on;
plot(mouthCoord(1), mouthCoord(2), 'x');
plot(ple(1), ple(2), 'o');
plot(pre(1), pre(2), 'o');
subplot(1,2,2);
imshow(eyeMapMouthFix);
hold on;
plot(mouthCoord(1), mouthCoord(2), 'x');


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

