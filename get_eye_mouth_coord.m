function [left_eye, right_eye, mouth] = get_eye_mouth_coord(img)
%Input image is the original RGB-image

eyePic = eyeMap(img);
mouthPic = mouthMap(img);

mouthPic = imadjust(mouthPic,stretchlim(mouthPic),[]);
mouthPic = uint8(rescale(mouthPic,0,255));

onlyEyes = uint8(face_threshold(eyePic, 200));
onlyMouth = uint8(face_threshold(mouthPic, 200)); 

logicalMouth = logical(onlyMouth);
logicalMouth = bwareafilt(logicalMouth,1);
s = regionprops(logicalMouth, 'centroid');
mouthCoord = s.Centroid;


logicalEyes = logical(onlyEyes);
[rows, cols] = size(eyePic);

%Zoom in on eye area using mouthcoordinates
logicalEyes( round(mouthCoord(2)) - round((1/8)*rows) : rows , : ) = 0; % bottom
logicalEyes( 1 : round((1/3.0)*rows), :) = 0;                  % top
logicalEyes( : , 1 : round(mouthCoord(1)) - round((1/5)*cols) ) = 0;    % left
logicalEyes( : , round(mouthCoord(1)) + round((1/4)*cols) : cols ) = 0; % right 
  
logicalEyes = bwareafilt(logicalEyes,2); %selects 2 largest objects
s = regionprops(logicalEyes, 'centroid');
eyeCoords = cat(1,s.Centroid); %Vector containg left and right eye position


left_eye = eyeCoords(1,:);
right_eye = eyeCoords(2,:);
mouth = mouthCoord;


