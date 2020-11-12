function [thresholded] = face_threshold(face, t_value)
%Thresholds the eyeMap
[rows, cols] = size(face);
thresholded = face;

for row = 1:rows
    for col = 1:cols
        if(face(row,col) > t_value)
            thresholded(row,col) = 255;
        else
            thresholded(row,col) = 0;
        end
    end
end
