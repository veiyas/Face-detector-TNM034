function [normalized_color_img] = whitePatch(rgb_img)
%Computes a white world normalization version of the rgb image
%Curently not giving any results 
redC = rgb_img(:,:,1);
greenC = rgb_img(:,:,2);
blueC = rgb_img(:,:,3);

rMax= max(redC, [], 'all'); %? sensor assumed as channel
gMax = max(greenC, [], 'all');
bMax = max(blueC, [], 'all');
%Seems kinda wacked, as it almost always returns 255 for each channel

alpha = gMax/rMax;
beta = gMax/bMax;

corrected_red = redC * alpha;
corrected_blue = blueC * beta;

normalized_color_img = cat(3, corrected_red, greenC, corrected_blue);

end

