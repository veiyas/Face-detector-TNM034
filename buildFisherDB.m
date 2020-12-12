function DB = buildFisherDB(dimensions)
c = 16; % Number of persons
%% Prepare training data
[images, numImages, correctIds] = loadImages('DB1_AND_DB2');

normalizedImages = cell(1, numImages);
for k = 1:numImages
    normalizedFace = normalizeFace(images{k});
    normalizedImages{k} = im2double(normalizedFace);
end

imgSize = size(normalizedImages{1});
imageVectors = zeros(imgSize(1) * imgSize(2), numImages);
for k = 1:numImages
    imageVectors(:,k) = normalizedImages{k}(:);
end

%% Calculate PCA to use later
globalMeanFace = mean(imageVectors, 2);
A = imageVectors - globalMeanFace;

[PCAEigVectors, PCAEigValues] = eig(A' * A, 'vector');

% Sort eigenvectors by eigenvalues
[~,ind] = sort(PCAEigValues, 'descend');
PCAEigVectors = PCAEigVectors(:,ind);

% Remove eigenvectors that are not needed
PCAEigVectors = PCAEigVectors(:,1:numImages-c);

% Get eigenvectors of the covariance matrix
Wpca = A * PCAEigVectors;
Wpca = Wpca ./ sqrt(sum(Wpca.^2,1));

PCAFaceCoords = transpose(Wpca) * A;

%% Find means and things on WPCA
vectorDim = size(PCAFaceCoords, 1);
numImagesOfPersons = zeros(c, 1);
meanOfPersons = zeros(vectorDim, c);
globalMean = zeros(vectorDim, 1);
% Each cell of groupedImages contains a matrix where the columns are
% different representations of the same person
groupedImages = cell(1, c);

for i = 1:numImages
    globalMean = globalMean + PCAFaceCoords(:,i);
    id =correctIds(i);

    numImagesOfPersons(id) = numImagesOfPersons(id) + 1;
    groupedImages{id}(:,numImagesOfPersons(id)) = PCAFaceCoords(:,i);
    
    meanOfPersons(:,id) = meanOfPersons(:,id) + PCAFaceCoords(:,i);
end

globalMean = globalMean ./ numImages;
meanOfPersons = meanOfPersons ./ numImagesOfPersons';

%% The LDA stuff

% Between-class scatter matrix
SB = zeros(vectorDim);
for i = 1:c
    diff = meanOfPersons(:,i) - globalMean;
    SB = SB + numImagesOfPersons(i) * (diff * diff');
end

% Within-class scatter matrix
SW = zeros(vectorDim);
for i = 1:c
   for j = 1:size(groupedImages{i}, 2)
       diff = groupedImages{i}(:,j) - meanOfPersons(:,i);
       SW = SW + (diff * diff');
   end
end


% Dubbelkolla detta
[Wfld, eigs] = eig(SB, SW);
[~, idx] = sort(diag(eigs), 'descend');
Wfld = Wfld(:, idx);
Wfld = Wfld(:, 1:dimensions);

Wopt = transpose(Wfld' * Wpca');

%% Find coordinates and save DB

% Should Wopt be normalized?

faceSpaceCoords = transpose(Wopt) * imageVectors;

% Create DB
DB.type = 'Fisher';
DB.meanFace = globalMeanFace;
DB.faceSpaceBasis = Wopt;
DB.faceSpaceCoords = faceSpaceCoords;
DB.ids = correctIds;

save DB;

%% Some random vis stuff that is very hard coded and might not work anymore?
% height = 400;
% width = 300;
% 
% imshow(reshape(globalMean, height, width)); title('Mean of all images');
% visualizeAll(meanOfPersons); title('Within-class means');
% 
% function visualizeAll(imgVec)
%     for iii = 1:size(imgVec, 2)
%         iiimages{iii} = reshape(imgVec(:,iii), height, width);
%         iiimages{iii} = rescale(iiimages{iii}, 0, 1);
%     end
%     montage(iiimages);
% end

end