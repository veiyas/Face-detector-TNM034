function [left_eye, right_eye, mouth] = get_eye_mouth_coord(img)
%Input image is the original RGB-image

eyePic = eyeMap(img);
YCBCR_im = rgb2ycbcr(img);
mouthPic = mouthMap(YCBCR_im);

J = imadjust(mouthPic,stretchlim(mouthPic),[]);
scaledJ = rescale(J,0,255);

T = graythresh(scaledJ);
T_uint8 = 255*T;

onlyEyes = face_threshold(eyePic, 200);
onlyMouth = face_threshold(scaledJ, 200); 

onlyEyes = uint8(onlyEyes);
onlyMouth = uint8(onlyMouth);

logicalMouth = logical(onlyMouth);
mouth = bwareafilt(logicalMouth,1); %selects 1 largest objects (magic)

mouthPos = regionprops( mouth, 'centroid');

mouthCoord = mouthPos.Centroid;

logicalEyes = logical(onlyEyes); %binary version of mouthMap
[rows, cols] = size(eyePic);

%CROPPING 
%Uncommented croppig code is "copied", change to imcrop later
 logicalEyes( round(mouthCoord(2)) - round((1/8)*rows) : rows , : ) = 0; % bottom
 logicalEyes( 1 : round((1/3.0)*rows), :) = 0;                  % top
 logicalEyes( : , 1 : round(mouthCoord(1)) - round((1/5)*cols) ) = 0;    % left
 logicalEyes( : , round(mouthCoord(1)) + round((1/4)*cols) : cols ) = 0; % right

% logicalEyes = imcrop(logicalEyes, []) %xmin, ymin, width, heigth
  
eyes = bwareafilt(logicalEyes,2); %selects 2 largest objects

eyePos = regionprops( eyes, 'centroid');
eyeCoords = cat(1,eyePos.Centroid); %Vector containg left and right eye position

left_eye = eyeCoords(1,:);
right_eye = eyeCoords(2,:);
mouth = mouthCoord;


