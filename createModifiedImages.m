function modifiedImages = createModifiedImages(image)
%GETMODIFIEDIMAGES Summary of this function goes here
%   ...

% "When testing your program for the examination, the facial images will be
% deliberately modified, to different levels of difficulty. The
% modifications will include translation of the face in the image, rotation
% (max +/- 5 degrees), scaling (max +/- 10%) and change of tone values
% in the image (max +/- 30%)." -- Course information

%% Generate rotated images
rotations = -5:1:5; % in degrees
rotatedImages = cell(1, length(rotations));
for i = 1:length(rotations)
    rotatedImages{i} = imrotate(image, rotations(i), 'bicubic', 'crop');
end

%% Generate scaled images
scales = 0.9:0.02:1.1;
scaledImages = cell(1, length(scales));
for i = 1:length(scales)
    scaledImages{i} = imresize(image, scales(i), 'bicubic');
end

% TODO Modify tone values

%%
%modifiedImages = [rotatedImages scaledImages];
modifiedImages = [rotatedImages]; % TODO Use all kinds of modifications

% TODO Mix different modifications

end

