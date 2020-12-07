function mouthMap = mouthMap(img)
%MOUTHMAP. Input image is RGB. Gives binary image!

imgYcbcr = rgb2ycbcr(img);
Cb = imgYcbcr(:,:,2);
Cr = imgYcbcr(:,:,3);

numPixels = size(imgYcbcr, 1) * size(imgYcbcr, 2);

crSquared = rescale(double(Cr).^2, 0, 1);
crOverCb = rescale(double(Cr) ./ double(Cb), 0, 1);
eta = 0.95 * (sum(crSquared, 'all')/numPixels) / (sum(crOverCb, 'all')/numPixels);

mouthMap = crSquared .* (crSquared - eta * crOverCb).^2;

% New stuff below
faceMask = getFaceMask(grayworldcorrection(img));
% These things work pretty good as long as the face mask is decent
% The face mask isn't decent for 4 images in DB2, otherwise it works well

% Strels to clean up binary image
closer = strel('disk', 5);
denoiser = strel('disk', 3);
vertLinesRemover = strel('line', 20, 0);

% Adjustments
mouthPic = mouthMap .* faceMask;
mouthPic = imadjust(mouthPic,stretchlim(mouthPic),[]);
mouthPic = uint8(rescale(mouthPic,0,255));

% Morph image reasonable binary image
logicalMouth = mouthPic > 240;
mouthMorphed = imopen(logicalMouth, vertLinesRemover);
mouthMorphed = imopen(mouthMorphed, denoiser);
mouthMorphed = imclose(mouthMorphed, closer);
mouthMap = mouthMorphed;
end

