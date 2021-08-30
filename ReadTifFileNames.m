function [tifFolder, tifNames] = ReadTifFileNames(tifFile)

tifFolder = fullfile(tifFile);
dirOutput = dir(fullfile(tifFolder, '*.tif')); % .tif .tiff .jpg
dirFileTif = sortObj(dirOutput);
tifNames = {dirFileTif.name}';
end