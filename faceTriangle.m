function[leftEye, rightEye, mouthPos] = faceTriangle(eyeMask, mouthMask)
%Get face triangle by looking at potential eyes and mouth and combining
%them.

eyeThresh = 220;
mouthThresh = 200;

mouthMask = imadjust(mouthMask,stretchlim(mouthMask),[]);
mouthMask = uint8(rescale(mouthMask,0,255));

eyeMask = uint8(face_threshold(eyeMask, eyeThresh));
eyeMask = logical(eyeMask);
mouthMask = uint8(face_threshold(mouthMask, mouthThresh));

%Get all potential eye candidates
eyeStats = regionprops('table',eyeMask, 'Centroid', 'BoundingBox');
%Store all centroids in a vector
eyeCentroids = cat(1,eyeStats.Centroid); 
eyesDetected = size(eyeCentroids);

%Get all potential mouth candidates
mouthStats = regionprops('table', mouthMask, 'centroid', 'BoundingBox');
mouthCentroids = cat(1, mouthStats.Centroid);
mouthsDetected = size(mouthCentroids);

%Get image size
[imgHeight, imgWidth] = size(eyeMask);

%Filter the mouth candidates
mouths = [];
k = 1;
for i=1:mouthsDetected
    %The mouth should be wider than its height
    if(mouthStats.BoundingBox(i,3) > mouthStats.BoundingBox(i,4))
        %The mouth should not be in the upper part of the image
        %Maybe even bottom half?
        if(mouthCentroids(i,2) > imgHeight/2)
            mouths(k,:) = mouthCentroids(i,:);
            k = k+1;
        end
    end
end

%Remove all eyes that does not fulfill certain attributes
ple = []; %Potential left eyes
l = 1; %Left eye counter
pre = []; %Potential right eyes
r = 1; %Right eye counter

combo = []; %Stores each possible eye+mouth combo

%This solution only works if we find ONLY 1 MOUTH.
for i=1:size(mouths,1)
    mouth = mouths(i,:);
    for j=1:eyesDetected
        eye = eyeCentroids(j,:);
        %Must be above mouth and within resonable distance. 200 is just empirical
        if(eye(2) > mouth(2) && abs(eye(2)-mouth(2)) < 200 )
            %if left side of mouth and within resonable distance.
            if(eye(1) < mouth(1) && abs(eye(1)-mouth(1)) < 75)
                ple(l,:) = eye;
                l = l+1;
            %if instead right side of mouth and within resonable distance.  
            elseif (eye(1) > mouth(1) && abs(eye(1)-mouth(1)) < 75)
                pre(r,:) = eye;
                r = r+1;
            end            
        end
    end
    
    %Choose the eye pair with the least y-coord difference
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
end

%Somehow choose the best potential combo and return eyes and mouth

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
















