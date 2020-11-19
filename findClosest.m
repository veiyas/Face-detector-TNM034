function [closestId, resNorm] = findClosest(normalizedImageVector, DB)
%FINDCLOSEST 
imgInFaceSpace = transpose(DB.faceSpaceBasis) * (normalizedImageVector - DB.meanFace);
error = imgInFaceSpace - DB.faceSpaceCoords;
[resNorm, closestId] = min(sqrt(sum(error.^2,1)));
end

