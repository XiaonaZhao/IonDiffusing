% function TaS2_batch(expName, tifPath, Mask, begin, rate, saveRoute, zone)
function TaS2_batch(expName, tifPath, Mask, begin, rate, saveRoute, zone, Fs)
% ITO

tic

% Fs = 100; % sampling rate
% All raw images
Value.tifFile = tifPath;
Value.tifDir = dir(fullfile(Value.tifFile, '*.tiff'));

% Potential and ScanRate
% r = rate/1000;
% r = -r/Fs;
% potential1 = (0 : r : -1.1)';
% end1 = potential1(end);
% potential2 = ((-2.2-end1-r) : (-r) : 0)';
% potential1 = (-0.6 : r : -1.1)';
% end1 = potential1(end);
% potential2 = ((-2.2 - end1 - r) : (-r) : -0.6)';
Value.potential = potentialLine(rate, Fs, -0.3, -0.8); % 2c
% Value.potential = [potential1' potential2']';% 1c


% Valid image sequence
Value.validDir = Value.tifDir(begin.frame:(begin.frame+length(Value.potential)));
Value.begin = begin;

% for n = 1:length(Value.maskNames)
Mask = imread(Mask);
mask = ~Mask;

if sum(mask(:)) == 0
    return
end

points = ReadTifMaskPoint(Value.tifFile, Value.validDir, mask);

col = size(points, 2);
curve = zeros(size(points));
parfor ii = 1:1:col
    % curve(:, ii) = lowp(points(:, ii), 1, 36, 0.1, 20, Fs); % SPR, 20;
%     curve(:, ii) = lowp(points(:, ii), 2, 12, 0.1, 20, Fs); % Bright Field, CV;
    curve(:, ii) = lowp(points(:, ii), 2, 12, 0.1, 20, 100); % Bright Field, CV, Fs = 20;
end
% curve = points; % imfilter
clear points

X = (1:1:(size(curve, 1)-1))';
dcurve = zeros(size(curve, 1)-1, size(curve, 2));

parfor ii = 1:col
    dcurve(:, ii) = diff(curve(:, ii));
end

Value.outside = -figSketch(dcurve);
outside = Value.outside;
img = figure('color','w');
hold on
for ii = 1:2
    plot(X, outside(:, ii), '.k')
end
xlabel('Frames'); ylabel('\DeltaIntensity''');
title([expName ' Na_2SO_4 ROI, ' num2str(rate) ' mV/s'])
hold off
figPath = [saveRoute '\' expName '_' num2str(zone) '_roi'];
saveas(img, figPath, 'fig')
% end

img2 = figure('color','w');
hold on
% for n = 1:length(Value.maskNames)
% outside = Value.outside;
for ii = 1:2
    plot(Value.potential, outside(:, ii), '.k')
end
% end
xlabel('Potential/V'); ylabel('\DeltaIntensity''');
title([expName ' \DeltaIntensity'' with Potential, Na_2SO_4, ' num2str(rate) ' mV/s'])
hold off
figPath2 = [saveRoute '\' expName '_intensityVSpotential_roi_' num2str(zone)];
saveas(img2, figPath2, 'fig')


% % ROIMEAN
tif0 = double(imread(fullfile(Value.tifFile, Value.tifDir(1).name)));
for ii = Value.begin.frame:(Value.begin.frame+length(Value.potential))
    tif  = (double(imread(fullfile(Value.tifFile, Value.tifDir(ii).name))) - tif0)./tif0;
%     tif  = (double(imread(fullfile(Value.tifFile, Value.tifDir(ii).name))) - tif0);
    Value.ROImean((ii-Value.begin.frame+1), 1) = ROImean(tif, mask); % TaS2
    Value.ROIsum((ii-Value.begin.frame+1), 1) = ROIsum(tif, mask); % TaS2
%     Value.ROImean((ii-Value.begin.frame+1), 1) = -ROImean(tif, mask); % TiS2
end

temp = lowp(Value.ROImean, 2, 12, 0.1, 20, Fs); % TaS2
% temp = lowp(Value.ROImean, 2, 12, 0.1, 20, 100); % Low sampleRate
% temp = lowp(Value.ROImean, 5, 22, 0.1, 20, 100); % TiS2 is not good
% enough
Value.dROImean = -diff(temp); clear temp
% Value.dROImean = -diff(Value.ROImean); % imfilter
img3 = figure('color','w');
plot(Value.potential, Value.dROImean, 'k')
xlabel('Potential/V'); ylabel('\DeltaIntensity''');
title([expName ' Averagered Intensity'' with Potential, Na_2SO_4, ' num2str(rate) ' mV/s'])
hold off
figPath3 = [saveRoute '\' expName '_AveragedintensityVSpotential_roi_' num2str(zone)];
saveas(img3, figPath3, 'fig')

temp = lowp(Value.ROIsum, 2, 12, 0.1, 20, Fs); % Low sampleRate
Value.dROIsum = diff(temp); clear temp
img4 = figure('color', 'w');
plot(Value.potential, Value.dROIsum, 'k')
xlabel('Potential/V'); ylabel('\DeltaIntensity''');
title([expName ' Total Intensity'' with Potential, Na_2SO_4, ' num2str(rate) ' mV/s'])
hold off
figPath4 = [saveRoute '\' expName '_TotalintensityVSpotential_roi_' num2str(zone)];
saveas(img4, figPath4, 'fig')


close all

% Save Value
cellpath = [saveRoute '\' expName '_' num2str(zone) '.mat']; 
save(cellpath, 'Value', '-v7.3');

toc

end