function[leftEyeCoord, rightEyeCoord, mouthCoord] = faceTriangle(img)
%Get face triangle by looking at potential eyes and mouth and combining
%them.

gwIM = grayworldcorrection(img);

% Get image size
[imgHeight, imgWidth] = size(gwIM);

% Get mouth map and mouth coords
mouthmap = mouthMap(img);
mouthCoord = getMouthCoord(mouthmap);

eyeMapBinMorphed = eyeMap(img);

% Remove all data below the mouth
eyeMapMouthFix = eyeMapBinMorphed;
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
l = 1; % Left eye counter
pre = []; % Potential right eyes
preArea = []; % Potenatial right eye areas
r = 1; % Right eye counter

heightDiffTolerance = 50;
areaDiffTolerance = 300;
normDiffTolerance = 10;

% Filter out candidates that are:
%   Not on correct side of mouth
%   Too far apart vertically
%   Not in similar size
%   Vector to mouth should be roughly same length
for i = 1:eyesDetected-1
    for j = i+1:eyesDetected
        tmpEyeArea1 = eyeAreas(i);
        tmpEyeArea2 = eyeAreas(j);
        
        tmpEyeMtrx = [eyeCentroids(i,:); eyeCentroids(j,:)];
        
        [~,leftEyeIndex] = min(tmpEyeMtrx, [], 1);
        [~,rightEyeIndex] = max(tmpEyeMtrx, [], 1);
        
        tmpLeftEye = tmpEyeMtrx(leftEyeIndex(1),:);
        tmpRightEye = tmpEyeMtrx(rightEyeIndex(1),:);
        
        leftEyeMouthNorm = norm(mouthCoord - tmpLeftEye);
        rightEyeMouthNorm = norm(mouthCoord - tmpRightEye);
        
        % Filters
        correctSides = tmpLeftEye(1) < mouthCoord(1) && mouthCoord(1) < tmpRightEye(1);
        belowMouth = tmpLeftEye(2) < mouthCoord(2) && tmpRightEye(2) < mouthCoord(2);
        smallVertDiff = abs(tmpLeftEye(1,2) - tmpRightEye(1,2)) < heightDiffTolerance;
        smallAreaDiff = abs(tmpEyeArea1 - tmpEyeArea2) < areaDiffTolerance;
        smallMouthNormDiff = abs(leftEyeMouthNorm - rightEyeMouthNorm) < normDiffTolerance;
        
        if correctSides && belowMouth && smallVertDiff && smallMouthNormDiff
            ple(l,:) = tmpLeftEye;
            pleAreas(l) = tmpEyeArea1;
            l = l+1;
            pre(r,:) = tmpRightEye;
            preAreas(r) = tmpEyeArea2;
            r = r+1;
        end
    end
end
% If we only have 1 candidate left its probably correct
% Else if choose the largest area ones
if size(ple,1) == 1 && size(pre,1) == 1
    leftEyeCoord = ple;
    rightEyeCoord = pre;
elseif size(ple,1) == 0 && size(pre,1) == 0
    % Guesstimated triangle, is this allowed?
    leftEyeCoord = mouthCoord - [100 100];
    rightEyeCoord = mouthCoord - [-100 100];
else 
    [~, maxAreaIndexLeft] = max(pleAreas);
    [~, maxAreaIndexRight] = max(preAreas);
    
    leftEyeCoord = ple(maxAreaIndexLeft,:);
    rightEyeCoord = pre(maxAreaIndexRight,:);
end
end








