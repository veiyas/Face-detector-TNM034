clear; clf;

% Should this be done in double or not?
input = im2double(imread('data/DB1/db1_01.jpg'));
inputYcbcr = rgb2ycbcr(input);
%imshow(inputYcbcr);

map = mouthMap(inputYcbcr);
imshow(map);