%% This program creates training and testing datasets based on CMEMS
%% 1. Create the original BayesNF training and testing datasets
load("SST_dataset.mat");
start_date = datetime('2016-1-1'); % Start date
end_date = datetime('2016-2-28');  % End date
date_indices = find(time_sst >= start_date & time_sst <= end_date);
time = time_sst(date_indices);
sla_dad = temp_cloud(:,:,date_indices);
temp_gt = temp_GT(:,:,date_indices);
aim_zos = sla_dad;

% Find positions that are NaN across all time (logical mask)
nan_mask = all(isnan(aim_zos), 3);   % true means all are NaN

% Assign these positions as 0
for k = 1:size(aim_zos,3)
    frame = aim_zos(:,:,k);
    frame(nan_mask) = 0;
    % frame(~nan_mask) = nan;
    aim_zos(:,:,k) = frame;
end
aim_lon = Lon_new;
aim_lat = Lat_new;
output_1322 = [ones(size(aim_zos,3),1)];
excel_lon = aim_lon;
excel_lat = aim_lat;
num_1322 = hours((time-datetime(1900,1,1))/24+2/24)';

% for i = 1:size(aim_zos)
%     subplot(1,2,1)
%     pcolor(aim_lon,aim_lat,aim_zos(:,:,i))
%     caxis([0,26])
%     shading flat
%     subplot(1,2,2)
%     pcolor(aim_lon,aim_lat,temp_gt(:,:,i))
%     caxis([0,26])
%     shading flat
%     pause
% end

count = 0;
for i = 1:size(aim_zos,1)
    for j = 1:size(aim_zos,2)
        count = count + 1;
        pro = [count.*ones(size(num_1322)),num_1322,excel_lon(i,j).*ones(size(num_1322)),excel_lat(i,j).*ones(size(num_1322)),squeeze(aim_zos(i,j,:))];
        output_1322 = [output_1322,pro];
    end
end
output_1322(:,1) = [];
excel_data_initial = [output_1322];
group = size(excel_data_initial,2)/5;
excel_data = [1,1,1,1,1];
for i = 1:group
    pro = excel_data_initial(:,1+5*(i-1):5*i);
    excel_data = [excel_data;pro];
end
excel_data(1,:) = [];

%% 2. Create Ground Truth
count = 0;
pro = [];
output_1322 = [ones(size(aim_zos,3),1)];
for i = 1:size(aim_zos,1)
    for j = 1:size(aim_zos,2)
        count = count + 1;
        pro = [count.*ones(size(num_1322)),num_1322,excel_lon(i,j).*ones(size(num_1322)),excel_lat(i,j).*ones(size(num_1322)),squeeze(temp_gt(i,j,:))];
        output_1322 = [output_1322,pro];
    end
end
output_1322(:,1) = [];
excel_data_initial = [output_1322];
group = size(excel_data_initial,2)/5;
excel_data_GT = [1,1,1,1,1];
for i = 1:group
    pro = excel_data_initial(:,1+5*(i-1):5*i);
    excel_data_GT = [excel_data_GT;pro];
end
excel_data_GT(1,:) = [];

%% 2. Remove NaN
location_nan = isnan(excel_data(:,5));

%% 1. Import prediction results
data_sequential = readtable('D:\Desktop\SST_train\output\SSIM_CMEMS_SST (2).csv');
pro_sla = table2array(data_sequential);

%% PICP
y = excel_data_GT(location_nan,5);
Q = pro_sla(location_nan,:);
% y = y(1:1,:);
% Q = Q(1:1,:);

%% Plotting
% Define confidence levels
alphas = [0.01:0.05:0.99];  % 1% to 99%
num_alpha = numel(alphas);

% Preallocate PICP
PICPs = zeros(1, num_alpha);

quantiles = 0.01:0.01:0.99;

% Loop through different confidence levels
for k = 1:num_alpha
    alpha = alphas(k);
    low_q  = (1-alpha)/2;
    high_q = 1 - (1-alpha)/2;
    
    [~, idx_low]  = min(abs(quantiles - low_q));
    [~, idx_high] = min(abs(quantiles - high_q));
    
    Qlow  = Q(:, idx_low);
    Qhigh = Q(:, idx_high);
    
    inside = (y >= Qlow) & (y <= Qhigh);
    PICPs(k) = mean(inside);
end

% Compute similarity metrics
MAE_curve  = mean(abs(PICPs - alphas));             % Mean Absolute Error
RMSE_curve = sqrt(mean((PICPs - alphas).^2));       % Root Mean Squared Error
ECE_mean   = mean(abs(PICPs - alphas));             % Expected Calibration Error

fprintf('Similarity between PICP and confidence levels:\n');
fprintf('MAE = %.4f, RMSE = %.4f, ECE = %.4f\n', MAE_curve, RMSE_curve, ECE_mean);

% Plot
figure('Position', [50, 50, 300, 400]);
plot(alphas*100, PICPs*100, 'r-o', 'LineWidth', 2);
hold on;
plot([0 100], [0 100], 'k--', 'LineWidth', 1.5);  % Ideal PICP reference line
xlabel('Confidence level α (%)');
ylabel('PICP (%)');
title('PICP vs Confidence Level');
legend('Inferred PICP', 'Ideal PICP', 'Location', 'Southeast');  % Legend
grid on;
xlim([0 100]);
ylim([0 100]);

% Annotate MAE, RMSE and ECE on the plot in percentage
text(5, 80, sprintf('MAE = %.2f%%', MAE_curve*100), 'FontSize', 11, 'Color', 'k');
text(5, 75, sprintf('RMSE = %.2f%%', RMSE_curve*100), 'FontSize', 11, 'Color', 'k');
text(5, 70, sprintf('ECE = %.4f', ECE_mean), 'FontSize', 11, 'Color', 'k');
