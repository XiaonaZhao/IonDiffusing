function OTM(tifPath, Mask)

tifDir = dir(fullfile(tifPath, '*.tiff'));
validDir = tifDir(begin.frame:(begin.frame + L -1));  % changed

tif0 = double(imread(fullfile(tifPath, validDir(1).name)));
tif = cell(L, 1); % images
tifValue = zeros(L, 1); % ROI value of the sequence
for ii = 1:L
    temp = double(imread(fullfile(tifPath, validDir(ii).name))) - tif0;
    temp = temp.*Mask;
    tif{ii, 1} = temp;
    tifValue(ii, 1) = sum(temp(:))/sum(Mask(:));
end

end