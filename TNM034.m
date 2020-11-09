%% normalizeFace testing DONT TOUCH
image = imread('data/DB1/db1_01.jpg');
image = im2double(image);

% Get eye and mouth points TEMPORARY SOLUTION
imshow(image);
hold on;
[xEye, yEye] = getpts;
[xMouth, yMouth] = getpts;
close;

% Plot face with triangle
faceTriangleX = [xEye' xMouth xEye(1)];
faceTriangleY = [yEye' yMouth yEye(1)];
subplot(1,2,1);
imshow(image);
hold on;
plot(faceTriangleX, faceTriangleY, 'red', 'LineWidth', 3);

% Rotate image by ensuring the line between the eyes is parallell to X-axis
eyeLineVec = [max(xEye) - min(xEye) max(yEye) - min(yEye) 0];
eyeLineVec = eyeLineVec ./ norm(eyeLineVec);

offsetAngle = atan2(norm(cross(eyeLineVec, [1 0 0])), dot(eyeLineVec, [1 0 0])) * (180 / pi);

subplot(1,2,2);
rotImage = imrotate(image, -offsetAngle, 'crop', 'bicubic');
imshow(rotImage)