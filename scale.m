function [scaledImage] = scale(mouthMapImage)

J = imadjust(mouthMapImage,stretchlim(mouthMapImage),[]);
scaledImage = rescale(J,0,255);

end

