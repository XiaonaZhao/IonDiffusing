function begin = triggerTime_2Cam(data, t, Fs_SPR, Fs_BF)
% check the timeline of double trigger

[~, begin.pikeSPR] = max(diff(data(1:30000, 1))); % trigger dot of the camera
[~, begin.pikeBF] = max(diff(data(1:30000, 4))); % trigger dot of the camera

plot(t(1:500000), data(1:500000, 2))
% plot(t, data(:, 2))

[x, ~] = ginput(2);
x1 = t((x(1)*10000):(x(2)*10000))';
y1 = data((x(1)*10000):(x(2)*10000), 2);
[x, ~] = ginput(2);
x2 = t((x(1)*10000):(x(2)*10000))';
y2 = data((x(1)*10000):(x(2)*10000), 2);

[px, ~] = crosspoint(x1, y1, x2, y2);
begin.CS1 = px*10000; % trigger dot of the CS studio


begin.beginSPR = ceil((begin.CS1 - begin.pikeSPR + 1)/10000*Fs_SPR);
begin.beginBF = ceil((begin.CS1 - begin.pikeBF + 1)/10000*Fs_BF);

end

