function TaS2(expName, maskPath, loader, rate, saveRoute, Fs)


[~, tifFile] = uigetfile('*.tiff', 'Multiselect', 'on', 'Read tif Folder');
tifDir = dir(fullfile(tifFile, '*.tiff'));

tic

if rate == 300
    r = -0.003;
    potential1 = (0 : r : -0.8)';
    potential2 = (-0.799 : (-r) : 0)';
    potential3 = (0.002 : r : -0.8)';
    potential4 = (-0.798 : (-r) : 0)';
    Value.potential = [potential1' potential2' potential3' potential4']';
    clear potential1 potential2 potential3 potential4
else
    
    r = -rate/(Fs*1000);
    
    potential1 = (0 : r : -0.8)';
    potential2 = ((-0.8-r) : (-r) : -0.0000001)';
    potential = [potential1' potential2' potential1' potential2']';
    clear potential1 potential2
end

varMat = load(loader);
begin = triggerTime(varMat.data, varMat.t);
validDir = tifDir(begin.frame:(begin.frame+length(potential)));

[~, maskNames] = ReadTifFileNames(maskPath);

hwait = waitbar(0, 'Please wait for the test >>>>>>>>');
for nn = 1:length(maskNames)
    Mask = imread(fullfile(maskPath, maskNames{nn}));
    mask = ~Mask;
    % span = 8;
    if sum(mask(:)) == 0
        return
    end
    
    points = ReadTifMaskPoint(tifFile, validDir, mask);
    
%     Fs = 100;
    col = size(points, 2);
    curve = zeros(size(points));
    for ii = 1:1:col
        curve(:, ii) = lowp(points(:, ii), 1, 36, 0.1, 20, Fs); % SPR, 20;
        %         curve(:, ii) = lowp(points(:, ii), 2, 11, 0.1, 20, Fs); % Bright Field, CV;
    end
    clear points
    
    X = (1:1:(size(curve, 1)-1))';
    dcurve = zeros(size(curve, 1)-1, size(curve, 2));
    
    for ii = 1:col
        dcurve(:, ii) = diff(curve(:, ii));
    end
    
    outside = figSketch(dcurve);
    img = figure('color','w');
    hold on
    for ii = 1:2
        plot(X, outside(:, ii), '.k')
    end
    xlabel('Frames'); ylabel('\DeltaIntensity''');
    title([expName ' Na_2SO_4 ROI' num2str(nn)])
    hold off
    figPath = [saveRoute '\' expName '_roi' num2str(nn) ];
    saveas(img, figPath, 'fig')
    
    processBar(length(Value.maskNames), nn, hwait)
    
end

cellpath = [saveRoute '\' expName '.mat'];
save(cellpath, 'Value', '-v7.3');

img2 = figure('color','w');
hold on
for nn = 1:length(Value.maskNames)
    outside = outside{nn, 1};
    for ii = 1:2
        plot(potential, outside(:, ii), '.k')
    end
end
xlabel('Potential/V'); ylabel('\DeltaIntensity''');
title([expName ' \DeltaIntensity'' with Potential, Na_2SO_4'])
hold off
figPath2 = [saveRoute '\' expName '_intensityVSpotential' num2str(nn) ];
saveas(img2, figPath2, 'fig')

toc

end