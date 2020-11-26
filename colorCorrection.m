function [normalized_color_img] = colorCorrection(rgb_img)
%Computes a grey world normalization version of the rgb image

redC = rgb_img(:,:,1);
greenC = rgb_img(:,:,2);
blueC = rgb_img(:,:,3);


red_avg = (sum(sum(redC)))./((size(redC,1)).*(size(redC,2)));
green_avg = (sum(sum(greenC)))./((size(greenC,1)).*(size(greenC,2)));
blue_avg = (sum(sum(blueC)))./((size(blueC,1)).*(size(blueC,2)));


alpha = green_avg/red_avg;
beta = green_avg/blue_avg;

corrected_red = redC * alpha;
corrected_blue = blueC * beta;

normalized_color_img = cat(3, corrected_red, greenC, corrected_blue);

end

