function begin = triggerTime(data, t, Fs)
% check the timeline of double trigger

[~, begin.pike] = max(diff(data(:, 1))); % trigger dot of the camera

plot(t, data(:, 2))
% plot(t, data(:, 2))

[x, ~] = ginput(2);
x1 = t((x(1)*10000):(x(2)*10000))';
y1 = data((x(1)*10000):(x(2)*10000), 2);
[x, ~] = ginput(2);
x2 = t((x(1)*10000):(x(2)*10000))';
y2 = data((x(1)*10000):(x(2)*10000), 2);

[px, ~] = crosspoint(x1, y1, x2, y2);
begin.CS1 = px*10000; % trigger dot of the CS studio

% % The following part is to check the accuracy of 'begin.CS1'
% [x, ~] = ginput(2);
% x3 = t((x(1)*10000):(x(2)*10000))';
% y3 = data((x(1)*10000):(x(2)*10000), 2);
% [x, ~] = ginput(2);
% x4 = t((x(1)*10000):(x(2)*10000))';
% y4 = data((x(1)*10000):(x(2)*10000), 2);
% 
% [px3, ~] = crosspoint(x3, y3, x4, y4);
% begin.CS3 = px3*10000; % recheck the trigger dot of the CS studio
% clear x1 x2 x3 x4 y1 y2 y3 y4 px px3

begin.frame = ceil((begin.CS1 - begin.pike + 1)/10000*Fs);

end

