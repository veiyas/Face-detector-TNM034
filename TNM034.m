%% normalizeFace ROTATION - DONT TOUCH
% Read and transform data
image = imread('data/DB1/db1_12.jpg');
image = im2double(image);
imsize = size(image);
imageYCBCR = rgb2ycbcr(image);

% Get eye and mouth points TEMPORARY SOLUTION
% This will be switched with eye + mouth maps
imshow(image);
[xEye, yEye] = getpts;
[xMouth, yMouth] = getpts;
close;

% Create triangle
faceTriangleX = [xEye' xMouth xEye(1)];
faceTriangleY = [yEye' yMouth yEye(1)];

% Get useful points
% The padded points are used to reduce interference from face features
leftmostEyePoint = min(faceTriangleX);
rightmostEyePoint = max(faceTriangleX);
bottommostMouthPoint = max(faceTriangleY);
topmostEyePoint = min(faceTriangleY);

% Rotation normalization
% Calculate eyes offset angle from x-axis
eyeLineVec = [max(xEye) - min(xEye) max(yEye) - min(yEye) 0];
eyeLineVec = eyeLineVec ./ norm(eyeLineVec);
offsetAngle = atan2(norm(cross(eyeLineVec, [1 0 0])), dot(eyeLineVec, [1 0 0]));

% Scaling normalization
% Get inverse CR channel with max contrast
CR = imageYCBCR(:,:,3);
stretchCR = imadjust(CR, stretchlim(CR), []);
invCR = 1 - stretchCR; 

% Do useful filtering
lplogFilter = fspecial('log',[15 15], 0.25);
lpLogGaussImage = imfilter(imgaussfilt(invCR, 2.0), lplogFilter);

% Remove data inside of face feature bounding box as it might interfere
lpLogGaussImage(topmostEyePoint-50:bottommostMouthPoint+20, leftmostEyePoint:rightmostEyePoint) = 0;

% Morhology to filter out unwanted details and accentuate the wanted edges
faceTopBotFinder = strel('line', 20, 0);
faceSideFinder = strel('line', 30, 90);
faceSides = imopen(lpLogGaussImage, faceSideFinder);
faceTopBot = imopen(lpLogGaussImage, faceTopBotFinder);

% Prepare morphed images
maxFaceSides = max(faceSides);
maxFaceTopBot = max(faceTopBot, [], 2)';

% Find left and right side of face
[~,faceLeftSide] = max(maxFaceSides(1:leftmostEyePoint));
[~,faceRightSide] = max(maxFaceSides(rightmostEyePoint:imsize(2)));
faceRightSide = faceRightSide + rightmostEyePoint; % Compensate for lost indices

% Find top and bot side of face
[~,faceTopSide] = max(maxFaceTopBot(1:topmostEyePoint));
[~,faceBotSide] = max(maxFaceTopBot(bottommostMouthPoint:length(maxFaceTopBot)));
faceBotSide = faceBotSide + bottommostMouthPoint; % Compensate for lost indices

% % Plot morphed images if needed
% subplot(1,3,1);
% imshow(lpLogGaussImage);
% subplot(1,3,2);
% imshow(faceSides);
% subplot(1,3,3);
% imshow(faceTopBot);

scaledImage = image(faceTopSide:faceBotSide, faceLeftSide:faceRightSide, :);
rotatedImage = imrotate(scaledImage, -offsetAngle * 180 / pi, 'crop', 'bicubic');

figure;
subplot(2,2,1);
imshow(image);
hold on;
plot(faceTriangleX, faceTriangleY, 'red', 'LineWidth', 3);
subplot(2,2,2);
imshow(rotatedImage);
subplot(2,2,3)
imshow(faceSides);
subplot(2,2,4);
imshow(faceTopBot);

