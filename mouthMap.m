function mouthMap = mouthMap(inputYcbcr)
%MOUTHMAP 

%numPixels = size(inputYcbcr, 1) * size(inputYcbcr, 2);

crSquared = rescale(inputYcbcr(:,:,3).^2, 0, 1);
crOverCb = rescale(inputYcbcr(:,:,3) ./ inputYcbcr(:,:,2), 0, 1);
eta = 0.95 * sum(crSquared, 'all') / sum(crOverCb, 'all');

mouthMap = crSquared .* (crSquared - eta * crOverCb).^2;

end

