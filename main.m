%% Build and save database
clear
clc
format compact

disp('Building DB');
[DB1Images, ~, ~] = loadImages('DB1');
dimensions = 7; % The number of dimensions in face space
DB = buildDB(DB1Images, dimensions);

save DB.mat DB

%% Test tnm034(im)
clc

[DB1Images, DB1NumImages, DB1CorrectIds] = loadImages('DB1');
[DB2Images, DB2NumImages, DB2CorrectIds] = loadImages('DB2');
[DB2ImagesNoBlur, DB2NumImagesNoBlur, DB2CorrectIdsNoBlur] = loadImages('DB2_EXCLUDING_BLUR');

disp('Results with unmodified training images (DB1)');
testWithNonmodifiedImages(DB1Images, DB1NumImages, DB1CorrectIds);
disp(' ');

% disp('Results with modified training images (DB1)');
% testWithModifiedImages(DB1Images, DB1NumImages, DB1CorrectIds);
% disp(' ');

disp('Results with unmodified DB2 images EXCLUDING BLURRY ONES');
testWithNonmodifiedImages(DB2ImagesNoBlur, DB2NumImagesNoBlur, DB2CorrectIdsNoBlur);
disp(' ');

% Test unknown faces
disp('Results with unknown faces (DB0)');
[DB0Images, DB0NumFiles] = loadImages('DB0');
numOfUnknownIdentifiedFaces = 0;
for k = 1:DB0NumFiles
    if tnm034(DB0Images{k}) ~= 0
        numOfUnknownIdentifiedFaces = numOfUnknownIdentifiedFaces + 1;
    end
end
fprintf('\tNumber of unknown faces identified:\t%i\n', numOfUnknownIdentifiedFaces);
disp(' ');
%% Local functions

function testWithNonmodifiedImages(images, numImages, correctIds)
    numCorrIds = 0;
    for i = 1:numImages
        if tnm034(images{i}) == correctIds(i)
            numCorrIds = numCorrIds + 1;
        else
            if tnm034(images{i}) == 0
                disp('hdjfklsjä');
            end
%             tnm034(images{i})
%             correctIds(i)
%             disp(' ');
        end
    end

    correctToTotalRatio = numCorrIds / numImages;
    fprintf('\tCorrect to total ratio:\t%f\n', correctToTotalRatio);
end

function testWithModifiedImages(images, numImages, correctIds)
    totNumImages = 0;
    numCorrIds = 0;
    for i = 1:numImages
        %i
        modifiedImages = createModifiedImages(images{i});
        numModImages = length(modifiedImages);
        totNumImages = totNumImages + numModImages;
        for k = 1:numModImages
            if tnm034(modifiedImages{k}) == correctIds(i) % correct
                numCorrIds = numCorrIds + 1;
            end
        end
    end

    correctToTotalRatio = numCorrIds / totNumImages;
    fprintf('\tCorrect to total ratio:\t%f\n', correctToTotalRatio);
end



















