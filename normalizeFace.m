function [doesFaceExist, normalizedImage] = normalizeFace(image)
%normalize Normalizes inputImage using eyecoords
%   Rotate, scale and tone to normalize

[leftEye, rightEye, mouth, numEyes] = get_eye_mouth_coord(image);
doesFaceExist = numEyes == 2;

if doesFaceExist == false
    normalizedImage = [];
    return
end

imageSize = [400, 300];
fixedLeftEye = [100 100];
fixedRightEye = [200 100];
fixedMouth = [150 300];
fixedPositions = [fixedLeftEye; fixedRightEye; fixedMouth];

transformation = fitgeotrans([leftEye; rightEye; mouth], fixedPositions, 'similarity');
normalizedImage = imwarp(image, transformation,'OutputView',imref2d(imageSize));

% Use grayscale image for now, might use something else later
normalizedImage = rgb2gray(normalizedImage);
end

