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
load DB.mat
%profile on

[DB1Images, DB1NumImages, DB1CorrectIds] = loadImages('DB1');
[DB2Images, DB2NumImages, DB2CorrectIds] = loadImages('DB2');

disp('Results with unmodified training images (DB1)');
testWithNonmodifiedImages(DB1Images, DB1NumImages, DB1CorrectIds, DB);
disp(' ');

disp('Results with modified training images (DB1)');
testWithModifiedImages(DB1Images, DB1NumImages, DB1CorrectIds, DB);
disp(' ');

disp('Results with unmodified DB2 images');
testWithNonmodifiedImages(DB2Images, DB2NumImages, DB2CorrectIds, DB);
disp(' ');

disp('Results with modified DB2 images');
testWithModifiedImages(DB2Images, DB2NumImages, DB2CorrectIds, DB);
disp(' ');

%% Test unknown faces
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

%profsave

%% Test specific parts of DB2
clear
load DB.mat
DBpart = ['DB2_bl'; 'DB2_cl'; 'DB2_ex'; 'DB2_il'];

disp(' ')
disp(' ======= PARTS OF DB2 ======= ')
disp(' ')

for i = 1:size(DBpart,1)
    disp(['Results for ' DBpart(i,:)]);
    [img, numImg, corrIds] = loadImages(DBpart(i,:));
    testWithNonmodifiedImages(img, numImg, corrIds, DB);
    disp(' ');
end

%% Local functions

function testWithNonmodifiedImages(images, numImages, correctIds, DB)
    numCorrIds = 0;
    falseRejection = 0;
    falseAccepted = 0;
    for i = 1:numImages
        classifiedAs = tnm034(images{i}, DB);
        if classifiedAs == correctIds(i)
            numCorrIds = numCorrIds + 1;
        else
            if classifiedAs == 0
                falseRejection = falseRejection + 1;
            else
                falseAccepted = falseAccepted + 1;
            end
        end
    end

    fprintf('\tCorrect / total:\t%i / %i\n', numCorrIds, numImages);
    fprintf('\tHow many faces were falesly rejected:\t%i\n', falseRejection);
    fprintf('\tHow many faces were falesly accepted:\t%i\n', falseAccepted);
    correctToTotalRatio = numCorrIds / numImages;
    fprintf('\tFRR:\t%f%%\n', (falseRejection/numImages)*100);
    fprintf('\tFAR:\t%f%%\n', (falseAccepted/numImages)*100);
    fprintf('\tAccuracy:\t%f%%\n', correctToTotalRatio*100);
end

function testWithModifiedImages(images, numImages, correctIds, DB)
    totNumImages = 0;
    numCorrIds = 0;
    falseRejection = 0;
    falseAccepted = 0;
    for i = 1:numImages
        %i
        modifiedImages = createModifiedImages(images{i});
        numModImages = length(modifiedImages);
        totNumImages = totNumImages + numModImages;
        for k = 1:numModImages
            classifiedAs = tnm034(modifiedImages{k}, DB);
            if classifiedAs == correctIds(i) % correct
                numCorrIds = numCorrIds + 1;
            else
                if classifiedAs == 0
                    falseRejection = falseRejection + 1;
                else
                    falseAccepted = falseAccepted + 1;
                end
            end
        end
    end

    fprintf('\tCorrect / total:\t%i / %i\n', numCorrIds, totNumImages);
    fprintf('\tHow many faces were falesly rejected:\t%i\n', falseRejection);
    fprintf('\tHow many faces were falesly accepted:\t%i\n', falseAccepted);
    correctToTotalRatio = numCorrIds / totNumImages;
    fprintf('\tFRR:\t%f%%\n', (falseRejection/totNumImages)*100);
    fprintf('\tFAR:\t%f%%\n', (falseAccepted/totNumImages)*100);
    fprintf('\tAccuracy:\t%f%%\n', correctToTotalRatio*100);
end



















