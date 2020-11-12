function [thresholded] = face_threshold(face)
%Thresholds the eyeMap
[rows, cols] = size(face);
thresholded = face;

for row = 1:rows
    for col = 1:cols
        if(face(row,col) > 210)
            thresholded(row,col) = 255;
        else
            thresholded(row,col) = 0;
        end
    end
end