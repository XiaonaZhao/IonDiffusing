prompt = {'Reference number of Data in Exp group:',...
    'Enter the folder of Masks for the Data:',...
    'Enter the matlab file of Data timeline:',...
    'Enter the scan rate of cyclic voltammetry (mV):',...
    'Save route:'};

dlg_title = 'Input';
num_lines = 1;

defaultans = {'A1',...
    'F:\TaS2\20190324_TaS2_0318_ITO\MaskA1',...
    'G:\TaS2\20190318_TaS2_ITO\matlab_data\B2.mat',...
    '100',...
    'F:\TaS2\20190324_TaS2_0318_ITO\result'};

answer = inputdlg(prompt, dlg_title, num_lines, defaultans);
expName = answer{1};
maskPath = answer{2};
loader = answer{3};
rate = str2num(answer{4});
saveRoute = answer{5};
 
% if DC
% TaS2_DC(expName, tifPath, maskPath, saveRoute);
% if CV
TaS2(expName, maskPath, loader, rate, saveRoute);

%%
load('G:\TaS2\TaS2_20190709_ITOx100\_Result\Wavelength_expTab.mat')

% for m = 1:size(expTab_TaS2, 1)
for m = [1 3 5 8]
    expName = expTab(m).expName;
    tifPath = expTab(m).tifPath;
    Mask = expTab(m).Mask;
    saveRoute = expTab(m).saveRoute;
    TaS2_DC(expName, tifPath, Mask, saveRoute);
end

%% dc
[~, A5.tifFile] = uigetfile('*.tif', 'Multiselect', 'on', 'Read tif Folder');
A5.tifDir = dir(fullfile(A5.tifFile, '*.tif'));
tif0 = double(imread(fullfile(A5.tifFile, A5.tifDir(1).name)));

highPotential = 0;
lowPotential = -1.1;
scanRate = 25;
sampleRate = 25;
potential = potentialLine(scanRate, sampleRate, highPotential, lowPotential);
A5.potential = potential;
A5.begin = begin.frame;

[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.*',  'All Files (*.*)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.roi',  'ROI (*.roi)'...
    },'Pick a .roi imageJ file');
[sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
A5.sROI = sROI;

A5.site = zeros(length(potential), length(sROI));
A5.location = cell(length(sROI), 1);


    
for ii = A5.begin:(A5.begin + length(potential) - 1)
    tif  = double(imread(fullfile(A5.tifFile, A5.tifDir(ii).name))) - tif0;
    
    for m = 1:length(sROI)
        [row, col] = ImageJroiLocation(sROI{m});
        img = tif((row(1):row(2)), (col(1):col(2)));
        A5.site((ii-A5.begin) + 1, m) = sum(img(:))/numel(img);
        A5.location{m, 1} = row;
    end
end

A5.cycle = (0:2/(length(potential)):2)';
A5.timeline = (1:length(potential))';

figure('color', 'w');
hold on
for ii = [1 length(sROI)]
% for ii = 1:length(sROI)
    temp = lowp(A5.site(:, ii), 2, 12, 0.1, 20, 100);
    % plot(A4.cycle(1:end-1), temp)
    plot(A5.timeline, temp)
end
xlabel('Frame');
ylabel('Intensity (a.u.)');
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off

A5.p1 = A5.location{1, 1};
A5.p2 = A5.location{end, 1};
A5.distance = norm(A5.p1-A5.p2); % Euclidean norm of vector

%%
ii = 26;
d = A5.part1{ii, 3};
d = imboxfilt(d, 3);
d(d > 0) = 0;
d = d.*mask;
imshow(d, 'DisplayRange',[], 'InitialMagnification', 'fit');
% colormap jet 
% impixelinfo

%%

A5.part1diff = cell(size(A5.part1, 1) -1, 1);
for ii = 1:length(A5.part1diff)
    A5.part1diff{ii} = A5.part1{ii+1, 2} - A5.part1{ii, 2};
end


temp = A5.part1diff;
A5.part1_3D = zeros([size(img) 349]);
for n = 1:length(temp)
    A5.part1_3D(:, :, n) = temp{n, 1};
end
clear temp

[A5.part1diff_value, A5.part1diff_time] = imgSeqDiff_s(A5.part1_3D);

%% 20190501_TaS2_ITO
figure('color', 'w');
hold on
for ii = 1:size(expTab, 1)
    plot(expTab(ii).num, expTab(ii).ROImean)
end
hold off

%% 20190501_TaS2_ITO
load('G:\TaS2\TaS2_1008_ITO\_Result\expTab.mat')

hwait = waitbar(0, 'Please wait for the test >>>>>>>>');
for m = 1:size(expTab, 1)
% for m = [3 4]
    expName = expTab(m).expName;
    tifPath = expTab(m).tifPath;
    Mask = expTab(m).Mask;
    begin = expTab(m).begin; 
    saveRoute = expTab(m).saveRoute;
    rate = expTab(m).ScanRate;
    zone  = expTab(m).zone;
%     Fs = 100;
    Fs = expTab(m).sampleRate;
    
    TaS2_batch(expName, tifPath, Mask, begin, rate, saveRoute, zone, Fs);
    
    disp(['The latest progress is about ' num2str(m) '.']);
    processBar(size(expTab, 1), m, hwait)
    
end
delete(hwait);

%%
figure('color', 'w');
hold on
plot(expTab(6).potential, expTab(6).dIntensity, 'color', [0.6350, 0.0780, 0.1840]) % 50 mV/s
plot(expTab(1).potential, expTab(1).dIntensity, 'color', [0.8500, 0.3250, 0.0980]) % 100 mV/s
plot(expTab(2).potential, expTab(2).dIntensity, 'color', [0.9290, 0.6940, 0.1250]) % 200 mV/s
plot(expTab(3).potential, expTab(3).dIntensity, 'color', [0.4660, 0.6740, 0.1880]) % 300 mV/s
plot(expTab(4).potential, expTab(4).dIntensity, 'color', [0.3010, 0.7450, 0.9330]) % 400 mV/s
plot(expTab(5).potential, expTab(5).dIntensity, 'color', [0, 0.4470, 0.7410]) % 500 mV/s
xlabel('Potential (V)');
ylabel('\DeltaIntensity'' (a.u.)');
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off
%%
figure('color', 'w');
hold on
plot(expTab(6).potential(3394:end,1), expTab(6).dROImean(3394:end,1), 'color', [0.6350, 0.0780, 0.1840]) % 50 mV/s
plot(expTab(1).potential(1698:end,1), expTab(1).dROImean(1698:end,1), 'color', [0.8500, 0.3250, 0.0980]) % 100 mV/s
plot(expTab(2).potential(850:end,1), expTab(2).dROImean(850:end,1), 'color', [0.9290, 0.6940, 0.1250]) % 200 mV/s
plot(expTab(3).potential(567:end,1), expTab(3).dROImean(567:end,1), 'color', [0.4660, 0.6740, 0.1880]) % 300 mV/s
plot(expTab(4).potential(426:end,1), expTab(4).dROImean(426:end,1), 'color', [0.3010, 0.7450, 0.9330]) % 400 mV/s
plot(expTab(5).potential(341:end,1), expTab(5).dROImean(341:end,1), 'color', [0, 0.4470, 0.7410]) % 500 mV/s
xlabel('Potential (V)');
ylabel('Averaged \DeltaIntensity'' (a.u.)');
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off
%%
figure('color', 'w');
hold on
plot(expTab(6).potential(3394:end,1), expTab(6).ROImean(3395:end,1), 'color', [0.6350, 0.0780, 0.1840]) % 50 mV/s
plot(expTab(1).potential(1698:end,1), expTab(1).ROImean(1699:end,1), 'color', [0.8500, 0.3250, 0.0980]) % 100 mV/s
plot(expTab(2).potential(850:end,1), expTab(2).ROImean(851:end,1), 'color', [0.9290, 0.6940, 0.1250]) % 200 mV/s
plot(expTab(3).potential(567:end,1), expTab(3).ROImean(568:end,1), 'color', [0.4660, 0.6740, 0.1880]) % 300 mV/s
plot(expTab(4).potential(426:end,1), expTab(4).ROImean(427:end,1), 'color', [0.3010, 0.7450, 0.9330]) % 400 mV/s
plot(expTab(5).potential(341:end,1), expTab(5).ROImean(342:end,1), 'color', [0, 0.4470, 0.7410]) % 500 mV/s
xlabel('Potential (V)');
ylabel('Intensity (a.u.)');
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off
%%
figure('color', 'w');
hold on
plot(expTab(6).timeline, -expTab(6).ROImean, 'color', [0.6350, 0.0780, 0.1840]) % 50 mV/s
plot(expTab(1).timeline, -expTab(1).ROImean, 'color', [0.8500, 0.3250, 0.0980]) % 100 mV/s
plot(expTab(2).timeline, -expTab(2).ROImean, 'color', [0.9290, 0.6940, 0.1250]) % 200 mV/s
plot(expTab(3).timeline, -expTab(3).ROImean, 'color', [0.4660, 0.6740, 0.1880]) % 300 mV/s
plot(expTab(4).timeline, -expTab(4).ROImean, 'color', [0.3010, 0.7450, 0.9330]) % 400 mV/s
plot(expTab(5).timeline, -expTab(5).ROImean, 'color', [0, 0.4470, 0.7410]) % 500 mV/s
% [0, 0, 1]
xlabel('Frame');
ylabel('Intensity (a.u.)');
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off
%%
figure('color', 'w');
hold on
plot(expTab(6).cycle, -expTab(6).ROImean, 'color', [0.6350, 0.0780, 0.1840]) % 50 mV/s
plot(expTab(1).cycle, -expTab(1).ROImean, 'color', [0.8500, 0.3250, 0.0980]) % 100 mV/s
plot(expTab(2).cycle, -expTab(2).ROImean, 'color', [0.9290, 0.6940, 0.1250]) % 200 mV/s
plot(expTab(3).cycle, -expTab(3).ROImean, 'color', [0.4660, 0.6740, 0.1880]) % 300 mV/s
plot(expTab(4).cycle, -expTab(4).ROImean, 'color', [0.3010, 0.7450, 0.9330]) % 400 mV/s
plot(expTab(5).cycle, -expTab(5).ROImean, 'color', [0, 0.4470, 0.7410]) % 500 mV/s
xlabel('Potential (V)');
ylabel('Intensity (a.u.)');
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off
%%
figure('color', 'w');
hold on
plot(expTab(6).potential(3201:end,1), expTab(6).ROImean(3202:end,1), 'color', [0.6350, 0.0780, 0.1840]) % 50 mV/s
plot(expTab(1).potential(1601:end,1), expTab(1).ROImean(1602:end,1), 'color', [0.8500, 0.3250, 0.0980]) % 100 mV/s
plot(expTab(2).potential(801:end,1), expTab(2).ROImean(802:end,1), 'color', [0.9290, 0.6940, 0.1250]) % 200 mV/s
plot(expTab(3).potential(535:end,1), expTab(3).ROImean(536:end,1), 'color', [0.4660, 0.6740, 0.1880]) % 300 mV/s
plot(expTab(4).potential(401:end,1), expTab(4).ROImean(402:end,1), 'color', [0.3010, 0.7450, 0.9330]) % 400 mV/s
plot(expTab(5).potential(321:end,1), expTab(5).ROImean(322:end,1), 'color', [0, 0.4470, 0.7410]) % 500 mV/s
xlabel('Potential (V)');
ylabel('Intensity (a.u.)');
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off
%%
figure('color', 'w');
hold on
plot(expTab(6).num, expTab(6).ROImean, 'color', [0.6350, 0.0780, 0.1840]) % 50 mV/s
plot(expTab(1).num, expTab(1).ROImean, 'color', [0.8500, 0.3250, 0.0980]) % 100 mV/s
plot(expTab(2).num, expTab(2).ROImean, 'color', [0.9290, 0.6940, 0.1250]) % 200 mV/s
plot(expTab(3).num, expTab(3).ROImean, 'color', [0.4660, 0.6740, 0.1880]) % 300 mV/s
plot(expTab(4).num, expTab(4).ROImean, 'color', [0.3010, 0.7450, 0.9330]) % 400 mV/s
plot(expTab(5).num, expTab(5).ROImean, 'color', [0, 0.4470, 0.7410]) % 500 mV/s
xlabel('Potential (V)');
ylabel('Intensity (a.u.)');
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off
%% TaS2_0806_Aux60
figure('color', 'w');
hold on
plot(expTab(1).potential(3522:end,1), expTab(1).dROImean(3522:end,1), 'color', [0.8500, 0.3250, 0.0980]) % 100 mV/s
plot(expTab(2).potential(4402:end,1), -6*expTab(2).dROImean(4402:end,1), 'color', [0.9290, 0.6940, 0.1250]) % 200 mV/s
xlabel('Potential (V)');
ylabel('Averaged \DeltaIntensity'' (a.u.)');
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off
%% Test the
% figure('color', 'w');
curve = lowp(Tpoint, 1, 6, 0.1, 20, 100); % Bright Field, CV;
plot(X, Tpoint,  X, curve)
% set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
% xlim([6800 7800])
% ylim([100 450])
%%
exp = cell(8, 8);
fields = {'expName', 'zone', 'tifPath', 'Mask', 'begin', 'ScanRate', 'sampleRate', 'saveRoute'};
expTab = cell2struct(exp, fields, 2);

%%
prefix = ('E:\TaS2\20190501_TaS2_ITO\Timer\used\');
d = sortObj(dir([prefix, '*.mat']));
for ii = 1:size(expTab, 1)
    varMat = load([prefix, d(ii).name]);
    % varMat = load('G:\EDL\MoS2_20191212\_Timer\A1_data.mat');
    % Fs = expTab(8).sampleRate;
    Fs = 100;
    begin = triggerTime_AC(varMat.data, varMat.t, Fs);
%     begin = triggerTime_MoS2(varMat.data, varMat.t, Fs);
    % Fs_SPR = 50; Fs_BF = 40;
    % begin = triggerTime_2Cam(varMat.data, varMat.t, Fs_SPR, Fs_BF);
    % begin = triggerTime_DC(varMat.data, Fs);
    expTab(ii).begin = begin;
    expTab(ii).data = varMat.data;
end



%% TaS2_20190627_ITO
for ii = 1:size(expTab, 1)
    [val, num] = min(expTab(ii).dROImean);
    potential = expTab(ii).potential;
    expTab(ii).reVal = val;
    expTab(ii).rePoten = potential(num);
end
%% TaS2_20190627_ITO
for ii = 1:size(expTab, 1)
%     expTab(ii).ROImean_lowp = lowp(expTab(ii).ROImean, 2, 12, 0.1, 20, 100);
    expTab(ii).ROImean_drifted = driftBaseline(expTab(ii).timeline, expTab(ii).ROImean);
end

%% 20190627_TaS2_ITO
% load('E:\TaS2\TaS2_20190627_ITO\_Result\expTab_p.mat')
savePath = 'E:\TaS2\TaS2_20190627_ITO\_Result\PuzzlePics\';
%%
figure('color','white');

for ii = 1:size(expTab_p, 1)
    expName = expTab_p(ii).expName;
    tifPath = expTab_p(ii).tifPath;
    sROI = expTab_p(ii).sROI;
    
    tifDir = dir(fullfile(tifPath, '*.tiff'));
    tif0 = double(imread(fullfile(tifPath, tifDir(potentialframe(1, ii)).name)));
    
    for jj = 2:size(potentialframe, 1) % for tif0 is the first
        tif  = double(imread(fullfile(tifPath, tifDir(potentialframe(jj, ii)).name)));
        tif = tif - tif0;
        
        [row, col] = ImageJroiLocation(sROI);
        img = tif((row(1):row(2)), (col(1):col(2)));
        
        localSums = imboxfilt(img, 11);
        
        imshow(localSums, 'DisplayRange',[], 'InitialMagnification', 'fit');
        title([num2str(ii*2), ' mV/s ', num2str(jj-1)]);
        %         c = cool;
        %         c = JShine;
        %         c = SublimeVivid;
        c = parula;
        c = flipud(c);
        map = colormap(c);
        colorbar;
        set(gca, 'CLim', [-800 0]);
        h = colorbar;
        set(get(h, 'title'), 'string', 'Intensity (a.u.)', 'FontSize', 12);
        pause(0.1);
        saveas(gcf,[savePath, expName, '_', num2str(jj, '%02d'), '.tif']);
    end
    
end

%% for TPO ======== TaS2_ITO_0324 ==========
load('E:\TaS2\20190324_TaS2_0318_ITO\result\expTab_TaS2_0324.mat')
% [cstrFilenames, cstrPathname] = uigetfile(...
%     {'*.*',  'All Files (*.*)';...
%     '*.zip',  'Zip-files (*.zip)';...
%     '*.roi',  'ROI (*.roi)'...
%     },'Pick a .roi imageJ file');
% [sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
[row, col] = ImageJroiLocation(sROI{1});
tifFile = 'E:\TaS2\20190324_TaS2_0318_ITO\B1_Na_-0-5V_110mVpp_0-2Hz_60-40_20s_Pike106fps';
tifDir = dir(fullfile(tifFile, '*.tiff'));
tif0 = double(imread(fullfile(tifFile, tifDir(236).name)));
tif1 = double(imread(fullfile(tifFile, tifDir(237+1).name)));
tif2 = double(imread(fullfile(tifFile, tifDir(261+1).name))); % primer = 392;
tif3 = double(imread(fullfile(tifFile, tifDir(392+1).name)));
% tif1 = 237; tif2 = 246; tif3 = 2814;
tif1 = (tif1-tif0)./tif0;
tif2 = (tif2-tif0)./tif0;
tif3 = (tif3-tif0)./tif0;
img1 = tif1((row(1):row(2)), (col(1):col(2)));
img2 = tif2((row(1):row(2)), (col(1):col(2)));
img3 = tif3((row(1):row(2)), (col(1):col(2)));

% figure('color', 'w');
localSums = imboxfilt(img2, 7);
imshow(localSums, 'DisplayRange',[], 'InitialMagnification', 'fit');
c = parula;
c = flipud(c);
map = colormap(c);
colorbar;
set(gca, 'CLim', [-0.1 0]);
h = colorbar;
set(get(h, 'title'), 'string', 'Intensity (a.u.)', 'FontSize', 12);


%%
cc = jet(10);
figure('color', 'w');
hold on
% for ii = 1:10
for ii = [2 4 6 8 10]
% ii = 2;
%     expTab(ii).redCurrent1 = 1/((expTab(ii).hole1site - expTab(ii).hole1Begin)/20);
%     expTab(ii).redCurrent = expTab(ii).hole1value/((expTab(ii).hole1site - expTab(ii).hole1Begin)/20);

%     plot(expTab(ii).cycle(1:end-1), expTab(ii).dROImean, 'color', cc(ii,:))
%     plot(expTab(ii).timeline, -expTab(ii).ROImean_drifted, 'color', cc(ii,:))
%     plot(expTab(ii).potential(1602:3202), expTab(ii).dROImean(1602:3202), 'color', cc(ii,:))
    
    y = lowp(expTab(ii).ROImean_lowp_drifted, 2, 12, 0.1, 20, 100);
    plot(expTab(ii).cycle(1:end), -smooth(y, 50)/10000, 'color', cc(ii,:))
%     plot(expTab(ii).cycle(1:end), -expTab(ii).ROImean_lowp_drifted, 'color', cc(ii,:))
%     ii = ii + 1;
% %     disp(['The ' num2str(ii) 'th progress is about ' num2str(expTab(ii).redCurrent1) '.']);
%     disp(['The next progress is about ' num2str(ii) '.']);
end
hold off
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
set(gca, 'linewidth', 1.5) 


%%
for ii = 1:10
    expTab(ii).hole1value = expTab(ii).ROImean_drifted(expTab(ii).hole1site);
    disp(['The ' num2str(ii) ' progress is about ' num2str(expTab(ii).hole1value) '.']);
end


%% 20190506_TaS2_ITO

prefix = ('E:\TaS2\TaS2_1008_ITO\_Result\Pics\fig2020\');
d = sortObj(dir([prefix, '*.mat']));
for ii = 1:size(expTab, 1)
% for ii = [3 4]
    load([prefix, d(ii).name]);
%     m = length(Value.validDir);
%     expTab(ii).line = (0:(-0.8/m):-0.8)';
%     expTab(ii).dIntensity = Value.outside(:, 1:2);
%     expTab(ii).potential = Value.potential;
%     expTab(ii).ROImean = Value.ROImean;
%     expTab(ii).dROImean = Value.dROImean;
    expTab(ii).tifPage = Value.tifPage;
    
%     expTab(ii).cycle = (0:2/(length(Value.dROImean)):2)';
%     expTab(ii).timeline = (1:length(Value.ROImean))';
end


%% 20191008_TaS2_ITO

prefix = ('G:\TaS2\TaS2_1008_ITO\_Result\_Mask\B3\');
d = sortObj(dir([prefix, '*.tif']));
mask = cell(5, 1);
mask{1} = ~imread([prefix, d(1).name]);
mask_img = mask{1};
for ii = 2:size(d, 1)
    mask{ii} = ~imread([prefix, d(ii).name]);
    mask_img = mask_img + mask{ii};
end
mask_img = ~mask_img;
figPath = [prefix '_B3_roi.tif'];
imwrite(mask_img, figPath)


%%
x = zeros(size(expTab, 1), 1);
y = zeros(size(expTab, 1), 1);
for ii = 1:size(expTab, 1)
    x(ii, 1) = expTab(ii).ScanRate;
    y(ii, 1) = expTab(ii).rePoten;
end
%%
Y = zeros(size(X));
Y1 = zeros(size(X, 1), 1);
Z = zeros(size(X, 1), 1);
for ii = 1:size(X, 1)
    Y(ii, 1) = expTab(ii).ROImean(X(ii, 1));
    Y(ii, 2) = expTab(ii).ROImean(X(ii, 2));
    Y1(ii, 1) = Y(ii, 2) - Y(ii, 1);
    Z(ii, 1) = (Y1(ii, 1))/(X(ii, 2) - X(ii, 1));
end

%%
c1 = polyfit(scanRate,oxiPeak,1); d1 = polyval(c1, scanRate);
c2 = polyfit(scanRate,rePeak,1); d2 = polyval(c2, scanRate);
% c3 = polyfit(X1,Y_ox2,1); d3 = polyval(c3, X1);

figure('color', 'w');
% plot(X1, Y_re, '.', X1, Y_ox1, '.', X1, Y_ox2, '.');
plot(scanRate,oxiPeak, '.', scanRate,rePeak, '.');
hold on
plot(scanRate, d1, 'color', [0.6350, 0.0780, 0.1840]);
plot(scanRate, d2, 'color', [0, 0.4470, 0.7410]);
% plot(X1, d3, 'color', [0.9290, 0.6940, 0.1250]);
hold off

%%
f = fittype('a*t + b*t^0.5 +c','independent','t','coefficients',{'a','b','c'});
% c1 = fit(scanRate, oxiPeak, f);
% c2 = fit(scanRate, rePeak, f);
c3 = fit(scanRate, reCurrent, f);

figure('color', 'w');
plot(c3, scanRate, reCurrent);
% hold on
% plot(c2, scanRate, rePeak);
xlabel('Scan Rate (mV/s)')
ylabel('Current (V)')
xlim([1 10])
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2); 
set(gca, 'linewidth', 1.5)
% hold off

%%
figure('color', 'w');
hold on
plot(expTab(2).potential, expTab(2).dIntensity, 'color', [0, 0.4470, 0.7410]) % K
plot(expTab(4).potential, expTab(4).dIntensity, 'color', [0.3010, 0.7450, 0.9330]) % Na
plot(expTab(6).potential, expTab(6).dIntensity, 'color', [0.4660, 0.6740, 0.1880]) % Li
plot(expTab(8).potential, expTab(8).dIntensity, 'color', [0.9290, 0.6940, 0.1250]) % Mg
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off

figure('color', 'w');
hold on
plot(expTab(2).potential, expTab(2).dROImean, 'color', [0, 0.4470, 0.7410]) % K 3.31
plot(expTab(4).potential, expTab(4).dROImean, 'color', [0.3010, 0.7450, 0.9330]) % Na 3.58
plot(expTab(6).potential, expTab(6).dROImean, 'color', [0.4660, 0.6740, 0.1880]) % Li 3.82
plot(expTab(8).potential, expTab(8).dROImean, 'color', [0.9290, 0.6940, 0.1250]) % Mg 4.28
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off

%% Rename image Seq
folder_name = uigetdir;
folder = dir(folder_name);

oldname = cell(length(folder)-2, 1);
for ii = 3:length(folder)
   oldname{ii-2} = folder(ii).name;
end

newname = cell(length(oldname), 1);
for ii = 2:length(oldname)
   a = oldname{ii};
   b = str2num(a(12:end-6));
   c = num2str(b, '%06d');
   newname{ii} = ['B4_', c, '.tiff'];
   movefile([folder_name '\' oldname{ii}], [folder_name '\' newname{ii}])
end

%%
saveRoute = 'F:\TaS2\20190506_TaS2_ITO\Result\Pic';
% V = [-0.72 -0.701 -0.682 -0.64 -0.604 -0.575]';
V = [-0.686 -0.675 -0.631]';
n = zeros(size(expTab, 1), length(V));
m = zeros(size(expTab, 1), length(V));
for ii = [2 3 6]
    tifDir = dir(fullfile(expTab(ii).tifPath, '*.tiff'));
    tif0 = double(imread(fullfile(expTab(ii).tifPath, tifDir(1).name)));

    a = expTab(ii).potential;
    for jj = 1:length(V)
        b = a(:) - V(jj);
        [~, n(ii, jj)] = min(abs(b));
        m(ii, jj) = n(ii, jj) + expTab(ii).begin.frame +length(expTab(ii).potential);
        tif = double(imread(fullfile(expTab(ii).tifPath, tifDir(m(ii, jj)).name))) - tif0;
        tif1 = tif((row(1):row(2)), (col(1):col(2)));
        figPath = [saveRoute '\B' num2str(ii) '_' num2str(jj) '.tif'];
        imwrite(tif1, figPath)
    end
end

%%
shift = zeros(6, 2);
for ii = 1:6
    shift(ii, 1) = (charge{ii, 1}(1,2) - charge{ii, 1}(2,2))./(charge{ii, 1}(1,1) - charge{ii, 1}(2,1));
    shift(ii, 2) = (charge{ii, 1}(3,2) - charge{ii, 1}(4,2))./(charge{ii, 1}(3,1) - charge{ii, 1}(4,1));
end
%% 20200525
% load('G:\TaS2\20190501_TaS2_ITO\_Result_100_std\B_group.mat')
% load('G:\TaS2\TaS2_20190627_ITO\_Result\expTab_B.mat')
load('G:\TaS2\TaS2_1008_ITO\_Result\expTab.mat')

hwait = waitbar(0, 'Please wait for the test >>>>>>>>');

% CurrentCurve = cell(size(expTab, 1), 3);
% for m = 1:size(expTab, 1)
CurrentCurve = cell(8, 3);
for m = 1:8
    expName = expTab(m).expName;
    tifPath = expTab(m).tifPath;
    Mask = expTab(m).Mask;
    begin = expTab(m).begin; 
    saveRoute = expTab(m).saveRoute;
    rate = expTab(m).ScanRate;
    zone  = expTab(m).zone;
    Fs = 100; % high scan rate group
%     Fs = 20; % low scan rate group
    
    [potential, tifPage, curve] = TaS2_batch2020(expName, tifPath, Mask, begin, rate, saveRoute, zone, Fs);
    CurrentCurve{m, 1} = potential;
    CurrentCurve{m, 2} = tifPage;
    CurrentCurve{m, 3} = curve;
    
    disp(['The latest progress is about ' num2str(m) '.']);
    processBar(size(expTab, 1), m, hwait)
    
end
delete(hwait);


%% for 20190627_ITO
cc = jet(10);
% m = 0;
figure('color', 'w');
hold on 
% for ii = [6 1 2 3 4 5]
for ii = 1:10
%     m = m+1;
    potential = CurrentCurve{ii, 1};
    curve = CurrentCurve{ii, 3};
%     plot(potential(size(potential,1)/2 : end), curve(size(potential,1)/2 : end), 'color', cc(m,:))
%     plot(potential(1:size(potential,1)/2), curve(1:size(potential,1)/2), 'color', cc(m,:))
    plot(potential, curve, 'color', cc(ii,:))
end

% m = 0;
title([expName ' smoothed curve, Na_2SO_4 at different concentration'])
xlabel('Potential (V)'); 
ylabel('Current (a.u.)');
% set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
% set(gca, 'linewidth', 1.5)

% figure('color', 'w');
% hold on 
% for ii = [6 1 2 3 4 5]
%     m = m+1;
%     potential = CurrentCurve{ii, 1};
%     tifPage = CurrentCurve{ii, 2};
%     plot(potential(size(potential,1)/2 : end), smooth(tifPage(size(potential,1)/2 : end),3), 'color', cc(m,:))
% end
% title([expName ' Optics2Electrics, Na_2SO_4 at different concentration'])
% xlabel('Potential (V)'); 
% ylabel('Current (a.u.)');
% set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
% set(gca, 'linewidth', 1.5)

%%
cc = jet(10);
figure('color', 'white');
hold on
% for ii = 1:size(Current, 1)
for ii = [2 4 6 8 10]
    m = length(Current(ii).potential);
    x = (1:ceil(m/2))';
    y = Current(ii).Current(1:ceil(m/2));
%     y = lowp(y, 2, 12, 0.1, 20, 100);
    jj = Current(ii).Span;
    yy1 = Current(ii).alpha * smooth(x, y, jj, 'loess');
%     plot(Current(ii).potential(1: ceil(m/2)), yy1(1: ceil(m/2)), 'color', cc(ii, :))
    % plot(Current(ii).potential(ceil(m/2):end), yy1(ceil(m/2):end), 'color', cc(ii, :))
    
    x = (ceil(m/2)+1:m-10)';
    y = Current(ii).Current(ceil(m/2)+1:m-10);
    jj = Current(ii).span;
    yy2 = smooth(x, y, jj, 'loess');
%     plot(Current(ii).potential(ceil(m/2)+1:end-10), yy2, 'color', cc(ii, :))
    
    yy = [yy1; yy2];
    plot(Current(ii).potential(1:m-10), yy, 'color', cc(ii, :))
    
    Current(ii).Current_smoothed = yy;
    Current(ii).re = min(yy);
    Current(ii).ox = max(yy);
%     title([num2str(ii), ' smooth span ', num2str(jj)]);
end
xlim([-1.1 -0.6])
xlabel('Potential (V vs. Ag wire)')
ylabel('Current density (a.u.)')
legend
hold off


%%
cc = jet(10);
figure('color', 'white');
hold on
% for ii = 1:size(Current, 1)
for ii = [6 1 2 3 4 5]
    m = length(Current(ii).potential);
    x = (1:ceil(m/2))';
%     y = Current(ii).Current(1:ceil(m/2));
     y0 = lowp(expTab(ii).ROImean_lowp_drifted, 2, 12, 0.1, 20, 100);
     y = -smooth(y0, 50)/1000;
     y = y(1:ceil(m/2));
%     y = lowp(y, 2, 12, 0.1, 20, 100);
    jj = Current(ii).Span;
    yy1 = Current(ii).alpha * smooth(x, y, jj, 'loess');
%     plot(Current(ii).potential(1: ceil(m/2)), yy1(1: ceil(m/2)), 'color', cc(ii, :))
    % plot(Current(ii).potential(ceil(m/2):end), yy1(ceil(m/2):end), 'color', cc(ii, :))
    
    x = (ceil(m/2)+1:m-10)';
    y = Current(ii).Current(ceil(m/2)+1:m-10);
    jj = Current(ii).span;
    yy2 = smooth(x, y, jj, 'loess');
%     plot(Current(ii).potential(ceil(m/2)+1:end-10), yy2, 'color', cc(ii, :))
    
    yy = [yy1; yy2];
    plot(Current(ii).potential(1:m-10), yy, 'color', cc(ii, :))
    
    Current(ii).Current_smoothed = yy;
%     title([num2str(ii), ' smooth span ', num2str(jj)]);
end
xlim([-1.1 -0.6])
xlabel('Potential (V vs. Ag wire)')
ylabel('Current density (a.u.)')
legend
hold off

%% for 20190501
cc = parula_linear(6);
Scanrate = [2; 4; 6; 8; 10];
% Scanrate = [50; 100; 200; 300; 400; 500];
% concentration = [50; 100; 150; 200; 250; 300];
Ipeak_re = zeros(6, 1);
Ipeak_ox = zeros(6, 1);
figure('color', 'white');
hold on
nn = 0;
% ii = 1;
for ii = 1:6
    nn = nn + 1;
    m = length(expTab(ii).potential);
%     x = (ceil(m/2)-1 : m)';
    x = (ceil(m/2)+1 : m)'; % for concentration in 20191008
    y = lowp(expTab(ii).tifPage, 2, 36, 0.1, 20, 100);
%     y = -smooth(y0, 50);
%     y = y(ceil(m/2)-1 : m); % COncentration has selected the last cycle.
%     jj = 1;
    jj = expTab(ii).Span;
    yy1 = smooth(x, y, jj, 'loess');
    %     plot(Current(ii).potential(1: ceil(m/2)), yy1(1: ceil(m/2)), 'color', cc(ii, :))
    % plot(Current(ii).potential(ceil(m/2):end), yy1(ceil(m/2):end), 'color', cc(ii, :))
    
%     plot(expTab(ii).potential(ceil(m/2)-1 : m), yy1, 'color', cc(nn, :))
    plot(expTab(ii).potential(ceil(m/2)+1 : m), yy1, 'color', cc(nn, :))% COncentration has selected the last cycle.
    
    Ipeak_re(nn, 1) = min(yy1);
    Ipeak_ox(nn, 1) = max(yy1);
    
end
xlim([-0.8 0])
xlabel('Potential (V vs. Ag/AgCl)')
ylabel('Current density (a.u.)')
box on
legend
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
set(gca, 'linewidth', 2) 
hold off


%% =========== for 20190318 the WideField&CV ===========
load('G:\TaS2\20190318_TaS2_ITO\result\B2.mat')
% y = lowp(tifPage, 2, 11, 0.1, 20, 100);
% x = (y(1:1601)+y(1601:3201))/2;

y = Value.outside{1, 1};
y1 = (y(:, 1)+y(:, 3))/2;
y1 = y1(2401:3200) - mean(y1(2950:3200, 1)).*ones(800, 1);
y2 = (y(:, 2)+y(:, 4))/2;
y2 = y2(1601:2400) - mean(y2(1601:2100, 1)).*ones(800, 1);
y = [y2; y1];

% B2_half = (B2(1:1601)+B2(1601:3201))/2;
% B2_filt = lowp(B2_half, 30*0.001, 8, 0.01, 20, 100);
% y = diff(B2_filt);
% plot(Value.potential(1:1600), y) 

figure('color', 'w');
plot(Value.potential(1:1600), y) 
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
set(gca, 'linewidth', 2, 'FontSize', 16)
% xlabel('Potential (V vs. Ag/AgCl)'); 
% ylabel('Current density (a.u.)');

hold on
yyaxis right
plot(potential, 1e6*current) % from the excel file of CHI A9.
% G:\TaS2\_CHI\CHI-e_1T_TaS2_250mMPBNa\A9_0V -0-8V_100mVpS.xlsx
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
set(gca, 'linewidth', 2, 'FontSize', 16)
% ylabel('Current (A)');

% legend on
hold off


%%
f_re = fit(Scanrate, Ipeak_re, 'poly1');
f_ox = fit(Scanrate, Ipeak_ox, 'poly1');
% f_re = fit(concentration, Ipeak_re, 'poly1');
% f_ox = fit(concentration, Ipeak_ox, 'poly1');

figure('color', 'white');
hold on
plot(f_ox, Scanrate, Ipeak_ox, '^');
plot(f_re, Scanrate, Ipeak_re, 'v');
% plot(f_ox, concentration, Ipeak_ox, '^');
% plot(f_re, concentration, Ipeak_re, 'v');
box on
xlabel('Scan rate (mV/s)')
ylabel('I_p (a.u.)')
legend
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
set(gca, 'linewidth', 2)
hold off


%% ====================TaS2_20190520_ITO_AFM==============================
tif0 = double(imread(fullfile(tifFile, tifDir(1).name)));
% L = size(tifDir, 1);
% L = 500;
% X = (1:L)';
% tif = cell(L, 1);
%% set the mask
I1 = tif0;
I2 = im2bw(mat2gray(I1), 261*0.001);
imshow(I2);
I2 = ~I2;
%%
% G:\TaS2\TaS2_20190520_ITO_AFM\Result\_sROI\RoiSetD2.zip
[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.*',  'All Files (*.*)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.roi',  'ROI (*.roi)'...
    },'Pick a .roi imageJ file');
[sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
[row1, col1] = ImageJroiLocation(sROI{1});
[row2, col2] = ImageJroiLocation(sROI{2});


%% view all data
% evr = 5;
% L = length(tifDir);
L = 400;
stage = zeros(L, 3); % 1, Delta Intensity; 2, Count of white dots; 3, Sum of white dots
fitcurve = cell(L, 2); % 1, histfit, fit x2; 2, histfit, fit y2;
for ii = 2:L
    % ii = 800;
    
    I3 = double(imread(fullfile(tifFile, tifDir(ii).name)));
    % I4 = double(imread(fullfile(tifFile, tifDir(ii+1).name)));
    % I = (imboxfilt(I3, evr) - tif0).*I2;
    I = round((tif0 - I3).*I2); % all
    % I = I(row1(1):row1(2), col1(1):col1(2)); % D2, below;
    % I = I(row2(1):row2(2), col2(1):col2(2)); % D2, above;
    % I = round(I2.*imboxfilt(I3 - tif0, evr));
    % temp = I3 - tif0;
    % temp = I4 - I3;
    % temp = round(imboxfilt((I4 - I3), 3));
    % I = temp;
    b = unique(I(:));
    b = b(b > 0);
    data = I(I < 0 | I > 0);
    % I = int8(temp);
    % imhist(I);
    % histogram(I, b(b > 0)); % histogram(C,Categories) plots only the subset of categories specified by Categories.
    % histogram(data, b)
    % hf = histfit(data, max(b), 'kernel');
    hf = histfit(data, length(b), 'kernel'); % 'kernel', 'normal'
    xlim([-800 800])
    % histogram(I, length(b));
    title(num2str(ii));
    
    % the data of the fitting curve
    x2 = get(hf(2),'XData');
    y2 = get(hf(2),'YData');
    
    [~, position] = max(y2);
    stage(ii, 1) = x2(position);
    stage(ii, 2) = length(data);
    stage(ii, 3) = sum(data);
    
    fitcurve{ii, 1} = x2;
    fitcurve{ii, 2} = y2;
    
    % histogram(I(:), 1000);
    % xlim([0 500])
    % ylim([0 2000])
    pause(0.01);
end

%%
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2); 
set(gca, 'linewidth', 2, 'FontSize', 16)

%%
c = fire(L);
line = (1:L)';
t = cell(L, 1);
for ii = 1:L
    t{ii, 1} = line(ii).*ones(100, 1)/100;
end

hold on
for ii = 2:L
    % plot3(fitcurve{ii, 1}, fitcurve{ii, 2}, t{ii, 1}, 'color', c(ii, :)); % [0.3010, 0.7450, 0.9330]
    plot3(t{ii, 1}, fitcurve{ii, 1}, fitcurve{ii, 2}, 'color', c(ii, :)); % [0.3010, 0.7450, 0.9330]
    
end
title('kernel') % 'kernel', 'normal'
xlabel('time (s)')
ylabel('\DeltaIntensity')
zlabel('Count')
hold off

%%
% the data of the histogram
% x1 = get(hf(1),'XData'); 
% y1 = get(hf(1),'YData');
 
% the data of the fitting curve
x2 = get(hf(2),'XData'); 
y2 = get(hf(2),'YData');

[~, position] = max(y2);
stage = x2(position);

%%
ii = 152;
I3 = double(imread(fullfile(tifFile, tifDir(ii).name)));
I4 = double(imread(fullfile(tifFile, tifDir(ii-1).name)));
temp = (I3 - I4).*I2;
line = temp(:);
tag = unique(line);
count = zeros(length(tag), 1);
for jj = 1:length(tag)
    count(jj, 1) = length(find(line == tag(jj)));
end
plot(tag, count)
xlim([5 350])



%%
whitedot = zeros(L-1, 1);
evr = 11;

for ii = 1:L-1
    I3 = double(imread(fullfile(tifFile, tifDir(ii+1).name)));
    % I4 = double(imread(fullfile(tifFile, tifDir(ii).name)));
    tif{ii, 1} = (imboxfilt(I3, evr) - tif0).*I2;
    % tif{ii, 1} = (imboxfilt(I3, evr) - imboxfilt(I4, evr)).*I2;
    temp = (tif{ii, 1} <= -20);
    whitedot(ii, 1) = sum(temp(:));
end
pause(0.05);
% figure('color', 'w');
plot(X(2:end), whitedot)


%% TaS2_20190327_B1 data
[~, tifFile] = uigetfile('*.tiff', 'Multiselect', 'on', 'Read tif Folder');
tifDir = dir(fullfile(tifFile, '*.tiff'));
tif0 = double(imread(fullfile(tifFile, tifDir(1).name)));
[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.*',  'All Files (*.*)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.roi',  'ROI (*.roi)'...
    },'Pick a .roi imageJ file');
[sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
figure('color', 'w');
[row, col] = ImageJroiLocation(sROI{1, 1});

savepath = 'G:\TaS2\20190324_TaS2_0318_ITO\result\B1_v\';
% 238-255, 451-455
tif = cell(25, 1);
for ii = 238:255
    tif1 = double(imread(fullfile(tifFile, tifDir(ii).name)));
    tif1 = (tif1-tif0)./tif0;
    tif{ii-237, 1} = tif1(row(1):row(2), col(1):col(2));
end

for ii = 450:456
    tif1 = double(imread(fullfile(tifFile, tifDir(ii).name)));
    tif1 = (tif1-tif0)./tif0;
    tif{ii-431, 1} = tif1(row(1):row(2), col(1):col(2));
end

for ii = 1:size(tif, 1)
    localSums = imboxfilt(tif{ii, 1}, 5);
    localSums(localSums > 0) = 0;
    imshow(localSums, 'DisplayRange',[], 'InitialMagnification', 'fit');
    if ii < 4 || ii > 21
        title(['t = ',  num2str((ii-1)*0.01), ' s,Voltage = 0 V'], 'FontSize', 18);
    else
        title(['t = ',  num2str((ii-1)*0.01), ' s, Voltage = -0.8 V'], 'FontSize', 18);
    end
    c = parula; % turbo parula fire
    c = flipud(c);
    colormap(c);
    set(gca, 'CLim', [-0.1 0]);
    h = colorbar;
    set(h, 'FontSize', 16);
    % set(get(h,'title'),'string','Intensity (a.u.)', 'FontSize', 16);
    
    pause(0.01);
    saveas(gcf,[savepath, 'B1_' num2str(ii, '%02d'), '.tif']);
end

tifFolder = 'G:\TaS2\20190324_TaS2_0318_ITO\result\B1_v\';
[~, tifNames] = ReadTifFileNames(tifFolder); %======.tif or .tiff=======
video = VideoWriter('TaS2_20190327_B1.avi'); % Prepare a .avi file
video.FrameRate = 3; % The same as the fps
open(video);
for ii = 1:size(tifNames, 1)  % The frame number for video
    frame = imread(fullfile(tifFolder, tifNames{ii}));
    writeVideo(video, frame);
end
close(video);


%%
tifFolder = 'D:\Projects\TaS2 - ions\histimg\histimg\';
[~, tifNames] = ReadTifFileNames(tifFolder); %======.tif or .tiff=======
video = VideoWriter('TaS2_20190520_D2.avi'); % Prepare a .avi file
video.FrameRate = 24; % The same as the fps
open(video);
for ii = 1:3:size(tifNames, 1)  % The frame number for video
    frame = imread(fullfile(tifFolder, tifNames{ii}));
    writeVideo(video, frame);
end
close(video);