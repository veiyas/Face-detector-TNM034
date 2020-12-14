function[leftEyeCoord, rightEyeCoord, mouthCoord] = faceTriangle(img)
%Get face triangle by looking at potential eyes and mouth and combining
%them.
gwIM = grayworldcorrection(img);

% Get image size
[imgHeight, imgWidth] = size(gwIM);

% Get mouth map and mouth coords
mouthmap = mouthMap(gwIM);
mouthCoord = getMouthCoord(mouthmap);

eyeMapBinMorphed = eyeMap(gwIM);

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
    mouthCoord = [428 382];
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
    
    leftEyeCoord = ple(lowestAngleIndex,:);
    rightEyeCoord = pre(lowestAngleIndex,:);
end
%[leftEyeCoord, rightEyeCoord] = findPupils(img, leftEyeCoord, rightEyeCoord);
end








