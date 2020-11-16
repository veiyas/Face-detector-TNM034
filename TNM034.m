%% DEBUGGING: Check normalization of all DB1 images
scaledImage = zeros(401,301,16);
jpgString = '.jpg';
beginString = 'data/DB1/db1_';
picIndexString = '';
for i = 1:16   
   if i < 10
       picIndexString = ['0' int2str(i)];
   else
       picIndexString = int2str(i);
   end   
    pathString = [beginString picIndexString jpgString];
    subplot(4,4,i);
    imshow(normalizeFace(imread(pathString)));   
end