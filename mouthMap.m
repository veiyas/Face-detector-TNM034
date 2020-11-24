function mouthMap = mouthMap(img)
%MOUTHMAP. Input image is RGB

imgYcbcr = rgb2ycbcr(img);
Cb = imgYcbcr(:,:,2);
Cr = imgYcbcr(:,:,3);

%numPixels = size(inputYcbcr, 1) * size(inputYcbcr, 2);

crSquared = rescale(double(Cr).^2, 0, 1);
crOverCb = rescale(double(Cr) ./ double(Cb), 0, 1);
eta = 0.95 * sum(crSquared, 'all') / sum(crOverCb, 'all');

mouthMap = crSquared .* (crSquared - eta * crOverCb).^2;

end

