function id = tnm034(im)
%
% im: Image of unknown face, RGB-image in uint8 format in the
% range [0,255]
%
% id: The identity number (integer) of the identified person,
% i.e. ‘1’, ‘2’,…,‘16’ for the persons belonging to ‘db1’
% and ‘0’ for all other faces.
%
% Your program code.
%%%%%%%%%%%%%%%%%%%%%%%%%%

% threshold = 1e-10;
threshold = 10;

% DB.mat should be precomputed to avoid extra caclulations
load DB.mat

% Should no detected faces be handled?
normalizedImg = im2double(normalizeFace(im));
imageVector = normalizedImg(:);

[idOfClosest, residualNorm] = findClosest(imageVector, DB);

if residualNorm < threshold
    id = idOfClosest;
else
    id = 0; % No face is close enough
end

end

% TO BE (RE?)MOVED
% %% DEBUGGING: Check normalization of all DB1 images
% scaledImage = zeros(401,301,16);
% jpgString = '.jpg';
% beginString = 'data/DB1/db1_';
% picIndexString = '';
% for i = 1:16   
%    if i < 10
%        picIndexString = ['0' int2str(i)];
%    else
%        picIndexString = int2str(i);
%    end   
%     pathString = [beginString picIndexString jpgString];
%     subplot(4,4,i);
%     imshow(normalizeFace(imread(pathString)));   
% end