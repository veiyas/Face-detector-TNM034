[DB2Images, DB2NumImages, DB2CorrectIds] = loadImages('DB2');
[DB1Images, DB1NumImages, DB1CorrectIds] = loadImages('DB1');

%%
imageNumDB2 = 16;

normalizedDB2 = normalizeFace(DB2Images{imageNumDB2});
normalizedDB1 = normalizeFace(DB1Images{DB2CorrectIds(imageNumDB2)});

subplot(131)
imshow(normalizedDB2)
subplot(132)
imshow(normalizedDB1)
subplot(133)
imshow(normalizedDB2 + normalizedDB1, []);