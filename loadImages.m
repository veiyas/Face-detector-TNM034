function [images, numImages] = loadImages(dataSet)
%LOADIMAGES Returns a cell array containing the images
%   dataSet is a string

if strcmp(dataSet, 'DB1')
    jpegFiles = dir('data/DB1/*.jpg');
elseif strcmp(dataSet, 'DB0')
    jpegFiles = dir('data/DB0/*.jpg');
end


numImages = length(jpegFiles);
images = cell(1, numImages);
for k = 1:numImages 
    images{k} = imread([jpegFiles(k).folder '/' jpegFiles(k).name]);
end

end

