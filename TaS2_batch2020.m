% function TaS2_batch(expName, tifPath, Mask, begin, rate, saveRoute, zone)
function [potential, tifPage, curve] = TaS2_batch2020(expName, tifPath, Mask, begin, rate, saveRoute, zone, Fs)
% ITO

tic

% Fs = 100; % sampling rate
% All raw images
Value.tifFile = tifPath;
Value.tifDir = dir(fullfile(Value.tifFile, '*.tiff'));

% Potential and ScanRate
potential = potentialLine(rate, Fs, 0, -0.8); % 2c high scan rate
% potential = potentialLine(rate, Fs, -0.6, -1.1); % 1c low scan rate
Value.potential = potential;

% Valid image sequence
Value.validDir = Value.tifDir(begin.frame:(begin.frame+length(Value.potential)));
Value.begin = begin;

% for n = 1:length(Value.maskNames)
Mask = imread(Mask);
mask = ~Mask;

if sum(mask(:)) == 0
    return
end

% In SI, the unit of charge, the coulomb, is defined as the charge carried
% by one ampere during one second.

tifPage = zeros(length(Value.potential), 1);
tif0 = double(imread(fullfile(Value.tifFile, Value.validDir(1).name)));
for ii = 1:length(Value.potential)
    tif1 = double(imread(fullfile(Value.tifFile, Value.validDir(ii).name)));
    tif2 = double(imread(fullfile(Value.tifFile, Value.validDir(ii+1).name)));
    tif = (tif2 - tif1)./tif0.*mask;
    tifPage(ii, 1) = sum(tif(:));
end
Value.tifPage = tifPage;

curve = lowp(tifPage, 2, 12, 0.1, 20, 100); % Bright Field, CV, Fs = 20;

% View plotting
img = figure('color', 'w');
plot(Value.potential, curve) % without smoothing or filtering the tifPage is too terrible
title([expName ' Optics2Electrics, Na_2SO_4, ' num2str(rate) ' mV/s'])
xlabel('Potential (V)'); 
ylabel('\DeltaIntensity (a.u.)');
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
set(gca, 'linewidth', 1.5)
figPath = [saveRoute '\' expName '_O2C_roi' num2str(zone)];
saveas(img, figPath, 'fig')
pause(0.1)


% Save Value
cellpath = [saveRoute '\2020_' expName '_' num2str(zone) '.mat']; 
save(cellpath, 'Value', '-v7.3');

close all
toc

end