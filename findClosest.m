function [closestId, resNorm] = findClosest(normalizedImageVector, DB)
%FINDCLOSEST

% For eigen
%imgInFaceSpace = transpose(DB.faceSpaceBasis) * (normalizedImageVector - DB.meanFace);
% For fisher
imgInFaceSpace = transpose(DB.faceSpaceBasis) * (normalizedImageVector);

error = imgInFaceSpace - DB.faceSpaceCoords;
[resNorm, idx] = min(sqrt(sum(error.^2,1)));
closestId = DB.ids(idx);
end

