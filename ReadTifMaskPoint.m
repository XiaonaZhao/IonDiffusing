function points = ReadTifMaskPoint(tifFile, tifDir, mask)

% [~, tifFile] = uigetfile('*.tiff', 'Multiselect', 'on', 'Read tif Folder');
% tifDir = dir(fullfile(tifFile, '*.tiff'));

row = size(tifDir, 1);

[m, ii] = find(mask ~= 0);
loc = [m ii];
col = sum(mask(:));
points = zeros(row, col);

tif0 = double(imread(fullfile(tifFile, tifDir(1).name)));
for ii = 1: row
%     tif = (double(imread(fullfile(tifFile, tifDir(ii).name))) - tif0).*mask;
    tif = (double(imread(fullfile(tifFile, tifDir(ii).name))) - tif0)./tif0.*mask;
    
    for jj = 1:col
        points(ii, jj) = tif(loc(jj, 1), loc(jj, 2));
    end
    
end
end