function [left_eye, right_eye, mouth, numberOfEyes] = get_eye_mouth_coord(img, faceMask)
%Input image is the original RGB-image

eyeThresh = 220;
mouthThresh = 200;

% eyePic = eyeMap(img);
eyePic = double(eyeMap(img)) .* faceMask;
eyePic = uint8(eyePic);
mouthPic = mouthMap(img);

mouthPic = imadjust(mouthPic,stretchlim(mouthPic),[]);
mouthPic = uint8(rescale(mouthPic,0,255));

onlyEyes = uint8(face_threshold(eyePic, eyeThresh));
onlyMouth = uint8(face_threshold(mouthPic, mouthThresh));

logicalMouth = logical(onlyMouth);
logicalMouth = bwareafilt(logicalMouth,1);
s = regionprops(logicalMouth, 'centroid');
mouthCoord = s.Centroid;

logicalEyes = logical(onlyEyes);
[rows, cols] = size(eyePic);

% %Zoom in on eye area using mouthcoordinates
% logicalEyes( round(mouthCoord(2)) - round((1/8)*rows) : rows , : ) = 0; % bottom
% logicalEyes( 1 : round((1/3.0)*rows), :) = 0;                  % top
% logicalEyes( : , 1 : round(mouthCoord(1)) - round((1/5)*cols) ) = 0;    % left
% logicalEyes( : , round(mouthCoord(1)) + round((1/4)*cols) : cols ) = 0; % right

logicalEyes = bwareafilt(logicalEyes,2); %selects 2 largest objects
s = regionprops(logicalEyes, 'centroid');
eyeCoords = cat(1,s.Centroid); %Vector containg left and right eye position

[r1,c1] = size(eyeCoords);
[r2,c2] = size(mouthCoord);

%If 2 Eyes and a mouth is found
if(r1 == 2 && r2 == 1)
    numberOfEyes = 2;
    left_eye = eyeCoords(1,:);
    right_eye = eyeCoords(2,:);
    mouth = mouthCoord;
%If only one eye is found, dont return eye or mouth position
%and set number of eyes to 1.
else
    left_eye = [];
    right_eye = [];
    mouth = [];
    if(r1 == 1) % If only one eye is found
        numberOfEyes = 1;
    else            %If no eyes are found
        numberOfEyes = 0;
    end
end

end



