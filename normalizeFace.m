function normalizedImage = normalizeFace(image)
%normalize Normalizes inputImage using eyecoords
%   Rotate, scale and tone to normalize

[leftEye, rightEye, mouth] = get_eye_mouth_coord(image);

eyeLineVec = [rightEye - leftEye 0];
eyeLineVec = eyeLineVec ./ norm(eyeLineVec);
offsetAngle = atan2(norm(cross(eyeLineVec, [1 0 0])), dot(eyeLineVec, [1 0 0]));

leftEye = round(leftEye);
rightEye = round(rightEye);
mouth = round(mouth);

% Size of scaled image
scaleY = 400;
scaleX = 300;
% Find scale factor to normalize eye position
targetEyesDist = 140; % This generally works well
distBetweenEyes = rightEye(1) - leftEye(1);
scaleFactor = targetEyesDist / distBetweenEyes;
image = imresize(image, scaleFactor, 'bicubic');

% Rescale all values
leftEye = leftEye .* scaleFactor;
rightEye = rightEye .* scaleFactor;
mouth = mouth .* scaleFactor;

topmostEye = min(leftEye(2), rightEye(2)); % Highest eye point
distEyesMouth = round((mouth(2) - topmostEye));
marginX = round((scaleX - distBetweenEyes)/2);

marginYtotal = scaleY - distEyesMouth;
marginYbot = round(0.35*marginYtotal);
marginYtop = round(0.65*marginYtotal);

rotatedImage = imrotate(image, offsetAngle * 180 / pi, 'crop', 'bicubic');
normalizedImage = rotatedImage(topmostEye-marginYtop:mouth(2)+marginYbot, leftEye(1)-marginX:rightEye(1)+marginX, :);
normalizedImage = imresize(normalizedImage, [scaleY scaleX]);
% Use grayscale image for now, might use something else later
normalizedImage = rgb2gray(normalizedImage);
end

