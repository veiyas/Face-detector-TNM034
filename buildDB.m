function DB = buildDB(imageVectors, dimensions)
%BUILDDB Returns DB created from array of normalized images in vector form

if dimensions > size(imageVectors, 2)
    error('Dimensions > number of images');
end

meanFace = mean(imageVectors, 2);
faceDiff = imageVectors - meanFace;

[eigVectors, eigValues] = eig(transpose(faceDiff) * faceDiff, 'vector');

% Sort eigenvectors by eigenvalues
[~,ind] = sort(eigValues, 'descend');
eigVectors = eigVectors(:,ind);

% Discard higher dimensions
eigVectors = eigVectors(:,1:dimensions);

% Get eigenvectors of the covariance matrix
covEigVectors = faceDiff * eigVectors;

% Normalize eigenvectors for use as basis for face space
faceSpaceBasis = covEigVectors ./ sqrt(sum(covEigVectors.^2,1));

% Coordinates in face space for the different faces
faceSpaceCoords = transpose(faceSpaceBasis) * faceDiff;

% Create DB
DB.meanFace = meanFace;
DB.faceSpaceBasis = faceSpaceBasis;
DB.faceSpaceCoords = faceSpaceCoords;
end

