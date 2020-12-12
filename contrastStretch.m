function [stretched] = contrastStretch(inputImage)
maxVal = max(inputImage(:));
minVal = min(inputImage(:));

stretched = (inputImage - minVal) / (maxVal - minVal);
end

