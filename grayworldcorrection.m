function [gwIM] = grayworldcorrection(rgb)
% Calculate lighting corrected image (gray world assumption)
alpha = mean(rgb(:,:,2), 'all') / mean(rgb(:,:,1), 'all');
beta = mean(rgb(:,:,2), 'all') / mean(rgb(:,:,3), 'all');

gwIM(:,:,1) = alpha .* rgb(:,:,1);
gwIM(:,:,2) = rgb(:,:,2);
gwIM(:,:,3) = beta .* rgb(:,:,3);
end

