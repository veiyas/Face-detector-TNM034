function DB = buildDB(images, dimensions)
%BUILDDB Returns DB created from images

%% Prepare training data
numImages = length(images);
normalizedImages = cell(1, numImages);
for k = 1:numImages
    normalizedImages{k} = im2double(normalizeFace(images{k}));
end
imgSize = size(normalizedImages{1});
imageVectors = zeros(imgSize(1) * imgSize(2), numImages);
for k = 1:numImages
    imageVectors(:,k) = normalizedImages{k}(:);
end

%% Be helpful :)
if dimensions > size(imageVectors, 2)
    error('Dimensions > number of images');
end

%% Dimensionality reduction
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

