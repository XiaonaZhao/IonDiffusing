function TaS2_DC(expName, tifPath, Mask, saveRoute)

Value.tifFile = tifPath;
Value.tifDir = dir(fullfile(Value.tifFile, '*.tiff'));

% Value.maskPath = maskPath;
% [~, Value.maskNames] = ReadTifFileNames(Value.maskPath);

% for n = 1:length(Value.maskNames)
n = 1;
    Mask = imread(Mask);
    mask = ~~Mask;
%     mask = ~Mask;
    % span = 8;
    if sum(mask(:)) == 0
        return
    end
    
%     points = ReadTifMaskPoint(Value.tifFile, Value.tifDir, mask);
%     
%     Fs = 106;
%     col = size(points, 2);
%     curve = zeros(size(points));
%     for ii = 1:1:col
%         curve(:, ii) = lowp(points(:, ii), 1, 36, 0.1, 20, Fs);
%     end
%     clear points
%     
    X = (1:1:size(Value.tifDir, 1))';
%     
%     Value.outside{n, 1} = figSketch(curve);
%     outside = Value.outside{n, 1};
%     img = figure('color','w');
%     for ii = 2:2:4
%         plot(X, outside(:, ii), '.k')
%     end
%     xlabel('Frames'); ylabel('\DeltaIntensity');
%     title([expName ' Na_2SO_4 ROI' num2str(n)])
%     figPath = [saveRoute '\' expName '_roi' num2str(n) ];
%     saveas(img, figPath, 'fig')
    
    tif0 = double(imread(fullfile(Value.tifFile, Value.tifDir(1).name)));
    roimean = zeros(size(X));
    for ii = 1:size(Value.tifDir, 1)
%         roi =  (double(imread(fullfile(Value.tifFile, Value.tifDir(ii).name))) - tif0).*mask;
        roi =  (double(imread(fullfile(Value.tifFile, Value.tifDir(ii).name))) - tif0)./tif0;
        roimean(ii) = -ROImean(roi, mask);
    end
    Value.roimean(:, n) = roimean;
    
    img0 = figure('color','w');
    plot(X, Value.roimean(:, n), 'k')
    xlabel('Frames'); ylabel('\DeltaIntensity');
    title([expName ' Na_2SO_4 averaged ROI' num2str(n)])
    figPath = [saveRoute '\' expName '_average_roi' num2str(n) ];
    saveas(img0, figPath, 'fig')
    
%     X = (1:1:(size(Value.tifDir, 1)-1))';
% %     dcurve = zeros(size(Value.tifDir, 1)-1, size(Value.tifDir, 2));
%     dcurve = diff(roimean);
% %     for ii = 1:col
% %         dcurve(:, ii) = diff(curve(:, ii));
% %     end
%     
%     Value.DiffOutside{n, 1} = figSketch(dcurve);
%     outside = Value.DiffOutside{n, 1};
%     img1 = figure('color','w');
%     hold on
%     for ii = 1:2
%         plot(X, outside(:, ii), '.k')
%     end
%     xlabel('Frames'); ylabel('\DeltaIntensity''');
%     title([expName ' Na_2SO_4 derivative ROI' num2str(n)])
%     hold off
%     figPath = [saveRoute '\' expName '_diff_roi' num2str(n) ];
%     saveas(img1, figPath, 'fig')
    
% end

cellpath = [saveRoute '\' expName '.mat'];
save(cellpath, 'Value', '-v7.3');

end