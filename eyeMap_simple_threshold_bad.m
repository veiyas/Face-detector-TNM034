function eyeMap_simple_threshold_bad(inputRGB)

%This method really only works on pictures with low background cluster!!!

RGB = imread(inputRGB);

gray = rgb2gray(RGB);

%imshow(gray)

%histogram(gray);

eqGray = histeq(gray);
imshow(eqGray)

EyeMap = eqGray <= 25;

imshow(EyeMap)