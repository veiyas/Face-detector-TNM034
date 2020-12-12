%% Build and save database
clear
clc
format compact

% disp('Building DB');
% [DB1Images, ~, ~] = loadImages('DB1');
% dimensions = 7; % The number of dimensions in face space
% DB = buildDB(DB1Images, dimensions);
% 
% save DB.mat DB
rank = 15;
buildFisherDB(rank);

%% Test tnm034(im)
clc

[DB1Images, DB1NumImages, DB1CorrectIds] = loadImages('DB1');
[DB2Images, DB2NumImages, DB2CorrectIds] = loadImages('DB2');
[DB2ImagesNoBlur, DB2NumImagesNoBlur, DB2CorrectIdsNoBlur] = loadImages('DB2_EXCLUDING_BLUR');

disp('Results with unmodified training images (DB1)');
testWithNonmodifiedImages(DB1Images, DB1NumImages, DB1CorrectIds);
disp(' ');

disp('Results with modified training images (DB1)');
testWithModifiedImages(DB1Images, DB1NumImages, DB1CorrectIds);
disp(' ');

disp('Results with unmodified DB2 images');
testWithNonmodifiedImages(DB2Images, DB2NumImages, DB2CorrectIds);
disp(' ');

disp('Results with modified DB2 images');
testWithModifiedImages(DB2Images, DB2NumImages, DB2CorrectIds);
disp(' ');

% disp('Results with unmodified DB2 images EXCLUDING BLURRY ONES');
% testWithNonmodifiedImages(DB2ImagesNoBlur, DB2NumImagesNoBlur, DB2CorrectIdsNoBlur);
% disp(' ');

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

%% Test specific parts of DB2
clear
DBpart = ['DB2_bl'; 'DB2_cl'; 'DB2_ex'; 'DB2_il'];

disp(' ')
disp(' ======= PARTS OF DB2 ======= ')
disp(' ')

for i = 1:size(DBpart,1)
    disp(['Results for ' DBpart(i,:)]);
    [img, numImg, corrIds] = loadImages(DBpart(i,:));
    testWithNonmodifiedImages(img, numImg, corrIds);
    disp(' ');
end

%% Local functions

function testWithNonmodifiedImages(images, numImages, correctIds)
    numCorrIds = 0;
    numNoId = 0; % how many times no face is identified
    for i = 1:numImages
        if tnm034(images{i}) == correctIds(i)
            numCorrIds = numCorrIds + 1;
        else
            if tnm034(images{i}) == 0
                numNoId = numNoId + 1;
            end
%                 subplot(121)
%                 imshow(images{i});
%                 subplot(122)
%                 imshow(normalizeFace(images{i}));
%                 title(['Incorrectly identified as ' num2str(tnm034(images{i}))]); % yes this is terribly inefficient
%                 pause;
                
                
%             tnm034(images{i})
%             correctIds(i)
%             disp(' ');
        end
    end

    fprintf('\tCorrect / total:\t%i / %i\n', numCorrIds, numImages);
    fprintf('\tHow many faces was not recognized:\t%i\n', numNoId);
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



















