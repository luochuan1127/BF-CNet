%% 1. Quantile output calibration
% s = 0.7705;
data_info = importdata('D:\Desktop\SST_train\output\SSIM_CMEMS_SST (2).csv');
pro_sla = data_info.data;

%% 2. Import validation dataset
data_info = importdata('D:\Desktop\SST_train\input\validation_CMEMS_SST.csv');
data_info = data_info.data;
pro_lon = data_info(:,4);
pro_lat = data_info(:,5);
pro_sla_initial = data_info(:,6);
pro_time_initial = data_info(:,3);
tab_time = tabulate(pro_time_initial);
tab_time = tab_time(tab_time(:,2)~=0,1);

%% 3. Import cloud mask data
load('date_list.mat')
load('input_sla_data.mat')
load('input_cloud_data.mat')
for i = 1:size(cloud_dad,3)
    pro_data = aim_zos_dad(:,:,i);
    loc_nan = isnan(cloud_dad(:,:,i));
    pro_data(loc_nan) = nan;
    sla_dad(:,:,i) = pro_data;
end
for i = 1:size(cloud_son,3)
    pro_data = aim_zos_son(:,:,i);
    loc_nan = isnan(cloud_son(:,:,i));
    pro_data(loc_nan) = nan;
    sla_son(:,:,i) = pro_data;
end
time_pro = datetime(string(date_list),'InputFormat','yyyyMMdd');

% Define date range
start_date = datetime('2016-1-1'); % Start date
end_date = datetime('2016-2-28');  % End date

% Find indices within date range
date_indices = find(time_pro >= start_date & time_pro <= end_date);
time = time_pro(date_indices);
sla_dad = sla_dad(:,:,date_indices);
sla_son = sla_son(:,:,date_indices);
aim_zos = sla_dad;
aim_lon = dad_lon;
aim_lat = dad_lat;
output_1322 = [ones(size(aim_zos,3),1)];
excel_lon = aim_lon;
excel_lat = aim_lat;
num_1322 = hours((time-datetime(1900,1,1))/24+2/24)';
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
location_nan = isnan(excel_data(:,5));

%% 3. Create matrix dataset
% Replace cloud region data with 0
% pro_sla(location_nan,:) = 0;
% pro_sla_initial(location_nan,:) = 0;

% Fill predicted results with available original data
% pro_sla(~location_nan,:) = repmat(pro_sla_initial(~location_nan,1), 1, 31);

for i = 1:length(tab_time)
    location = find(pro_time_initial==tab_time(i));
    matrix_data(i).lon = pro_lon(location);
    matrix_data(i).lat = pro_lat(location);
    matrix_data(i).sla_initial = pro_sla_initial(location);
    matrix_data(i).sla = pro_sla(location,:);
    matrix_data(i).cloud = location_nan(location,:);
end

%% 4. Extract target sea region
for i = 1:length(matrix_data)
    pro_sla_gt = matrix_data(i).sla_initial;
    matrix_sla_gt(:,:,i) = reshape(pro_sla_gt,71,86); % pro_sla(:,5) represents SLA values of day 5
    pro_cloud = matrix_data(i).cloud;
    matrix_cloud(:,:,i) = reshape(pro_cloud,71,86);   % pro_sla(:,5) represents SLA values of day 5
end
numDays = size(pro_sla,2);
matrix_sla_bayesnf = zeros(71,86, numDays, length(matrix_data)); % Pre-allocate 4D matrix
for i = 1:length(matrix_data)
    for j = 1:numDays
        pro_sla_bayesnf = matrix_data(i).sla(:, j);
        matrix_sla_bayesnf(:, :, j, i) = reshape(pro_sla_bayesnf, 71,86);
    end
end

%% 5. SSIM
% Build quantile list
quantiles = 0.01:0.01:0.99;  % Closed interval, step size 0.01
numSamples = size(matrix_sla_gt, 3);    % 1578
numQuantiles = size(matrix_sla_bayesnf, 3); % 31
SSIM_values = zeros(numSamples, numQuantiles);

% Compute SSIM for all samples and quantiles
for i = 1:numSamples
    mask = ~matrix_cloud(:,:,i);  % Mask of the i-th sample
    if ~any(mask(:))
        % Skip if all zero
        continue
    end

    gt_layer = matrix_sla_gt(:,:,i);
    gt_masked = gt_layer;
    gt_masked(~mask) = 0;  % Outside mask set to 0

    for j = 1:numQuantiles
        pred_layer = matrix_sla_bayesnf(:,:,j,i);
        pred_masked = pred_layer;
        pred_masked(~mask) = 0;  % Outside mask set to 0
        
        pred_masked(isnan(pred_masked)) = 0;
        gt_masked(isnan(gt_masked)) = 0;

        SSIM_values(i,j) = ssim(pred_masked, gt_masked); % Compute SSIM
    end
end

mean_SSIM = round(mean(SSIM_values), 6);
[max_SSIM, max_idx] = max(mean_SSIM);
best_quantile = quantiles(max_idx);

% Display results
fprintf('Maximum average SSIM: %.6f\n', max_SSIM);
fprintf('Corresponding quantile: %.2f\n', best_quantile);
fprintf('Quantile index: %d\n', max_idx);
