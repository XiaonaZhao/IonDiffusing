% cv2_min = min(cv2')';
cv3_min = min(cv3')';
cv3_draft = zeros(size(cv3));
for ii = 1:length(cv3_min)
    cv3_draft(ii, :) = cv3(ii, :) - cv3_min(ii, 1)*ones(1, size(cv3, 2));
end

x = (1:22)';
xi = (1:0.025:22)';
cv3_draft_intp = zeros(length(xi), size(cv3_draft, 2));
for ii = 1:size(cv3_draft, 2)
    Y = cv3_draft(:, ii);
    cv3_draft_intp(:, ii) = interp1q(x,Y,xi);
end
%%
cv3_1cy = cv3_draft_intp;
%%
figure('color','white');
imshow(cv3_1cy, 'DisplayRange',[], 'InitialMagnification', 'fit');
colormap jet
% figure('color','white');
% imshow(log(cv3_1cy(:, 1:end)), 'DisplayRange',[], 'InitialMagnification', 'fit');
% colormap jet

%%
cv3_diff = zeros(size(cv3_1cy, 1)-1, size(cv3_1cy, 2));
for ii = 1:size(cv3_1cy, 1)-1
    cv3_diff(ii, :) = cv3_1cy(ii+1, :) - cv3_1cy(ii, :);
end
%%
figure('color','white');
imshow(cv3_diff(:, 1:end), 'DisplayRange',[], 'InitialMagnification', 'fit');
colormap jet

%%
figure('color','white');
ii = 1;
% line = smooth(cv2(ii, :), 10);
plot(X, cv2(ii,:));
% findpeaks(cv2(ii,:), X, 'MinPeakDistance',20)
% xlim([50 600])
%%
peakinfo{ii, 1} = cursor_info;
%%
cv2_1cy_std = zeros(size(cv2_1cy_trans));
for ii = 1:size(cv2_1cy_trans, 2)
    cv2_1cy_std(:, ii) = (cv2_1cy_trans(:, ii) - a*min(cv2_1cy_trans(:, ii)))./(max(cv2_1cy_trans(:, ii))-min(cv2_1cy_trans(:, ii)));
end