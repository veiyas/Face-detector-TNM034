function [images, numImages, correctIds] = loadImages(dataSet)
%LOADIMAGES Returns a cell array containing the images
%   dataSet is a string used as some kind of enum

if strcmp(dataSet, 'DB1')
    jpegFiles = dir('data/DB1/*.jpg');
elseif strcmp(dataSet, 'DB0')
    jpegFiles = dir('data/DB0/*.jpg');
elseif strcmp(dataSet, 'DB2')
    jpegFiles = dir('data/DB2/*.jpg');
    
% Parts of DB2 below
elseif strcmp(dataSet, 'DB2_bl')
    jpegFiles = dir('data/DB2/bl*.jpg');
    
elseif strcmp(dataSet, 'DB2_cl')
    jpegFiles = dir('data/DB2/cl*.jpg');
    
elseif strcmp(dataSet, 'DB2_ex')
    jpegFiles = dir('data/DB2/ex*.jpg');
    
elseif strcmp(dataSet, 'DB2_il')
    jpegFiles = dir('data/DB2/il*.jpg');
    
elseif strcmp(dataSet, 'DB2_EXCLUDING_BLUR')
    jpegFiles = dir('data/DB2_EXCLUDING_BLUR/*.jpg');
end


numImages = length(jpegFiles);
images = cell(1, numImages);
correctIds = zeros(numImages, 1);
for k = 1:numImages
    
    if ~strcmp(dataSet, 'DB0')
        correctIds(k) = str2double(extractBetween(jpegFiles(k).name, "_", ".jpg"));
    end
    
    images{k} = imread([jpegFiles(k).folder '/' jpegFiles(k).name]);
end

end

