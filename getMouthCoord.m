function [mouthCoord] = getMouthCoord(mouthmap)
%Get all potential mouth candidates
[imgHeight, imgWidth] = size(mouthmap);
mouthStats = regionprops('table', mouthmap, 'Centroid', 'BoundingBox', 'Area');
mouthCentroids = cat(1, mouthStats.Centroid);
mouthAreas = cat(1, mouthStats.Area);
mouthsDetected = size(mouthCentroids);

% The region with the largest area is 99% the mouth
% So we search for it with some extra filters as precaution
mouthCoord = [];
currArea = -1;
k = 1;
widthLengthRatioLimit = 1.8;
for i = 1:mouthsDetected
    % The mouth should be wider than its height
    if(mouthStats.BoundingBox(i,3) > mouthStats.BoundingBox(i,4))
        % The mouth should not be in the upper part of the image
        % Maybe even bottom half?
        if(mouthCentroids(i,2) > imgHeight/2)
            % The region needs to be a fair bit wider than tall
            if mouthStats.BoundingBox(i,3) / mouthStats.BoundingBox(i,4) > widthLengthRatioLimit
                % The region needs to have a larger area
                if (mouthAreas(i) > currArea)
                    currArea = mouthAreas(i);
                    mouthCoord = mouthCentroids(i,:);
                    k = k+1;
                end
            end
        end
    end
end

if currArea == -1 % If no mouth candidate was found we use placeholder
    disp('Could not find mouth');
    mouthCoord = [441.5 415.6]; 
end

end

