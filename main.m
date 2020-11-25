%% Build and save database
clear
format compact

% Prepare training data
jpegFiles = dir('data/DB1/*.jpg'); 
numFiles = length(jpegFiles);
images = cell(1, numFiles);
normalizedImages = cell(1, numFiles);
for k = 1:numFiles 
    images{k} = imread(['data/DB1/' jpegFiles(k).name]);
    normalizedImages{k} = im2double(normalizeFace(images{k}));
end
imgSize = size(normalizedImages{1});
imageVectors = zeros(imgSize(1) * imgSize(2), numFiles);
for k = 1:numFiles 
    imageVectors(:,k) = normalizedImages{k}(:);
end

% Calculate DB
dimensions = 15; % The number of dimensions in face space
DB = buildDB(imageVectors, dimensions);

save DB.mat DB

%% Test tnm034(im)

% For now the training data is used unmodified
% TODO apply rotations and so on...
knownFacesIdentificationResults = zeros(numFiles,1);
for i = 1:numFiles
    knownFacesIdentificationResults(i) = tnm034(images{i});
end

knownFacesCorrectResults = (1:numFiles)'; % Not the most robust...
numCorrectIds = nnz(...
    knownFacesIdentificationResults == knownFacesCorrectResults);

correctToTotalRatio = numCorrectIds / numFiles

% DOES NOT WORK RIGHT NOW, eyeCoords in get_eye_m... gets empty,
% not sure why atm
% Test unknown faces
% jpegFilesDB0 = dir('data/DB0/*.jpg'); 
% numFilesDB0 = length(jpegFilesDB0);
% numOfUnknownIdentifiedFaces = 0;
% for k = 3:numFilesDB0 
%     image = imread(['data/DB0/' jpegFilesDB0(k).name]);
%     if tnm034(image) ~= 0
%         numOfUnknownIdentifiedFaces = numOfUnknownIdentifiedFaces + 1;
%     end
% end
% 
% numOfUnknownIdentifiedFaces

