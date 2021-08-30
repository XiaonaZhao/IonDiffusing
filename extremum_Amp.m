function output = extremum_Amp(sample)

frame = size(sample, 3);
row = size(sample, 1);
col = size(sample, 2);
output = zeros(row, col);
for ii = 1:row
    parfor jj = 1:col
        temp = sample(ii, jj, :);
        temp = reshape(temp, frame, 1);
        Ma = max(temp(:));
        Mi = min(temp(:));
        output(ii, jj) = Ma - Mi;
    end
end
end