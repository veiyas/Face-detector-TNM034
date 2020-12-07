%% Build and save database
clear
format compact

[images, numImages] = loadImages('DB1');
dimensions = 7; % The number of dimensions in face space
DB = buildDB(images, dimensions);

save DB.mat DB

%% Test tnm034(im)
clc

disp('Results with unmodified training images (DB1)');
testWithNonmodifiedImages(images, numImages)
disp(' ');

disp('Results with modified training images (DB1)');
testWithModifiedImages(images, numImages);
disp(' ');


% DOES NOT WORK RIGHT NOW, eyeCoords in get_eye_m... gets empty,
% not sure why atm
% Test unknown faces
disp('Results with unknown faces (DB0)');
jpegFilesDB0 = dir('data/DB0/*.jpg'); 
numFilesDB0 = length(jpegFilesDB0);
numOfUnknownIdentifiedFaces = 0;
for k = 1:numFilesDB0
    image = imread(['data/DB0/' jpegFilesDB0(k).name]);
%     if k == 4 || k == 1 % OOOOOBBS DETTA ÄR BARA FÖR ATT UNDVIKA KRASCH
%         continue;
%     end
    
    if tnm034(image) ~= 0
        numOfUnknownIdentifiedFaces = numOfUnknownIdentifiedFaces + 1;
    end
end
fprintf('\tNumber of unknown faces identified:\t%i\n', numOfUnknownIdentifiedFaces);
disp(' ');
%% Local functions

% Dont necessarily work for anything other than DB1
function testWithNonmodifiedImages(images, numImages)
    knownFacesIdentificationResults = zeros(numImages,1);
    
    for i = 1:numImages
        knownFacesIdentificationResults(i) = tnm034(images{i});
    end

    knownFacesCorrectResults = (1:numImages)'; % Not the most robust...
    numCorrectIds = nnz(...
        knownFacesIdentificationResults == knownFacesCorrectResults);

    correctToTotalRatio = numCorrectIds / numImages;
    fprintf('\tCorrect to total ratio:\t%f\n', correctToTotalRatio);
end

% Dont necessarily work for anything other than DB1
function testWithModifiedImages(images, numImages)
    totNumImages = 0;
    numCorrIds = 0;
    for i = 1:numImages
        %i
        modifiedImages = createModifiedImages(images{i});
        numModImages = length(modifiedImages);
        totNumImages = totNumImages + numModImages;
        for k = 1:numModImages
            if tnm034(modifiedImages{k}) == i % correct
                numCorrIds = numCorrIds + 1;
            end
        end
    end

    correctToTotalRatio = numCorrIds / totNumImages;
    fprintf('\tCorrect to total ratio:\t%f\n', correctToTotalRatio);
end



















