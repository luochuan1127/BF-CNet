%% 1. Quantile output calibration
data_sequential = readtable('D:\Desktop\SST_train\output\SSIM_CMEMS_SST (2).csv');
pro_sla = table2array(data_sequential);

%% 2. Import validation dataset
data_info = importdata('D:\Desktop\SST_train\input\validation_CMEMS_SST.csv');
data_info = data_info.data;
pro_lon = data_info(:,4);
pro_lat = data_info(:,5);
pro_sla_initial = data_info(:,6);
pro_time_initial = data_info(:,3);
tab_time = tabulate(pro_time_initial);
tab_time = tab_time(tab_time(:,2)~=0,1);

%% 3. Create matrices
for i = 1:length(tab_time)
    location = find(pro_time_initial==tab_time(i));
    matrix_data(i).lon = pro_lon(location);
    matrix_data(i).lat = pro_lat(location);
    matrix_data(i).sla_initial = pro_sla_initial(location);
    matrix_data(i).sla = pro_sla(location,:);
end

%% 4. Extract target sea region
for i = 1:length(matrix_data)
    pro_sla_gt = matrix_data(i).sla_initial;
    matrix_sla_gt(:,:,i) = reshape(pro_sla_gt,71,86); % pro_sla(:,5) represents SLA values of day 5
end
numDays = size(pro_sla,2);
matrix_sla_bayesnf = zeros(71, 86, numDays, length(matrix_data)); % Pre-allocate 4D matrix
for i = 1:length(matrix_data)
    for j = 1:numDays
        pro_sla_bayesnf = matrix_data(i).sla(:, j);
        matrix_sla_bayesnf(:, :, j, i) = reshape(pro_sla_bayesnf, 71, 86);
    end
end

%% 5. Visualization (MAE)
load("SST_dataset.mat"); 
nan_mask = isnan(matrix_sla_gt); 
all_nan_mask = all(nan_mask, 3); 
lon_range = [107,121]; 
lat_range = [8,25]; 
start_date = datetime('2016-1-1'); % Start date 
end_date = datetime('2016-2-28');  % End date 
date_indices = find(time_sst >= start_date & time_sst <= end_date); 
time = time_sst(date_indices); 
sla_dad = temp_cloud(:,:,date_indices); 

% Identify land areas: positions that are NaN across all dates
land_mask = all(isnan(temp_cloud), 3); % true = land

% Calculate missing-data ratio per day (excluding land)
perc = zeros(1, size(sla_dad, 3)); 
for i = 1:size(sla_dad,3) 
    frame = sla_dad(:,:,i); 
    ocean_mask = ~land_mask; 
    total_ocean = sum(ocean_mask(:)); 
    nan_ocean = sum(isnan(frame(ocean_mask))); 
    perc(i) = nan_ocean / total_ocean; 
end

% Offset annotation position (top-left corner)
dx = (lon_range(2)-lon_range(1))*0.02; 
dy = (lat_range(2)-lat_range(1))*0.02;

% ---------- Select days 21–27 ----------
day_range = 21:27;   
num_days = numel(day_range);

% Create compact layout: 3 rows × num_days columns
t = tiledlayout(3, num_days, 'TileSpacing','compact', 'Padding','compact');

for j = 1:num_days
    i = day_range(j);  

    % ---------- Row 1: cloud ----------
    nexttile(j)
    pro_3 = sla_dad(:,:,i);   
    m_proj('mercator','lon',lon_range,'lat',lat_range);
    m_pcolor(Lon_new,Lat_new,pro_3); 
    m_gshhs('lc','patch',[0.5 0.5 0.5],'Edgecolor','k'); 
    m_grid('xtick',[],'ytick',[],'tickdir','in','linest',':','xaxisloc','top');
    colormap(flipud(nclCM(150,20)))
    cloud_pct = round(perc(i)*100,1);
    title({datestr(time(i),'yyyy-mm-dd'), ['Clouds: ', num2str(cloud_pct), '%']})
    m_text(lon_range(1)+dx, lat_range(2)-dy, sprintf('(a-%d)', j), ...
           'color','w','fontweight','bold','fontsize',10, ...
           'vertical','top','horizontal','left');

    % ---------- Row 2: Ground Truth ----------
    nexttile(j+num_days)
    pro_2 = squeeze(matrix_sla_gt(:,:,i));
    pro_2(all_nan_mask) = nan;
    cloud_mask = isnan(pro_3);
    pro_2 = pro_2';
    pro_2(~cloud_mask) = nan;  
    m_proj('mercator','lon',lon_range,'lat',lat_range);
    m_pcolor(Lon_new,Lat_new,pro_2);
    m_gshhs('lc','patch',[0.5 0.5 0.5],'Edgecolor','k');
    m_grid('xtick',[],'ytick',[],'tickdir','in','linest',':','xaxisloc','top');
    colormap(flipud(nclCM(150,20)))
    m_text(lon_range(1)+dx, lat_range(2)-dy, sprintf('(b-%d)', j), ...
           'color','w','fontweight','bold','fontsize',10, ...
           'vertical','top','horizontal','left');

    % ---------- Row 3: BF-CNet ----------
    nexttile(j+2*num_days)
    pro_1 = squeeze(matrix_sla_bayesnf(:,:,49,i));
    pro_1(all_nan_mask) = nan;
    diff_map = pro_1';
    diff_map(~cloud_mask) = nan;  
    m_proj('mercator','lon',lon_range,'lat',lat_range);
    m_pcolor(Lon_new,Lat_new,diff_map);
    m_gshhs('lc','patch',[0.5 0.5 0.5],'Edgecolor','k');
    m_grid('xtick',[],'ytick',[],'tickdir','in','linest',':','xaxisloc','top');
    colormap(flipud(nclCM(150,20)))
    m_text(lon_range(1)+dx, lat_range(2)-dy, sprintf('(c-%d)', j), ...
           'color','w','fontweight','bold','fontsize',10, ...
           'vertical','top','horizontal','left');
end

% Manually create colorbar
cb = colorbar('Position',[0.96 0.2 0.015 0.4]); 
cb.Label.String = 'SST/℃';
cbarrow

% ---------- Adjust window size ----------
set(gcf,'Position',[216,-223,1600,940]);

%% Calculate cloud-region metrics
% for i = 5:11
% Pre-allocate storage vectors
MAE_all  = zeros(1,58);
RMSE_all = zeros(1,58);
SSIM_all = zeros(1,58);
PSNR_all = zeros(1,58);

for i = 1:58
    % -------- Bayesnf --------
    pro_1 = squeeze(matrix_sla_bayesnf(:,:,49,i));
    pro_1(all_nan_mask) = nan;

    % -------- GT --------
    pro_2 = squeeze(matrix_sla_gt(:,:,i));
    pro_2(all_nan_mask) = nan;

    % -------- Cloud --------
    pro_3 = sla_dad(:,:,i);

    % --- Construct mask (valid only in cloud regions) ---
    mask = isnan(pro_3') & ~isnan(pro_1) & ~isnan(pro_2);

    % Data after applying mask
    v1 = pro_1(mask);
    v2 = pro_2(mask);

    % -------- MAE & RMSE --------
    mae_val = mean(abs(v1 - v2));
    rmse_val = sqrt(mean((v1 - v2).^2));

    % -------- PSNR --------
    L = max([v1; v2]) - min([v1; v2]);
    psnr_val = 10*log10(L^2 / mean((v1 - v2).^2));

    % -------- SSIM --------
    tmp1_masked = zeros(size(pro_1));
    tmp2_masked = zeros(size(pro_2));
    tmp1_masked(mask) = v1;
    tmp2_masked(mask) = v2;
    ssim_val = ssim(tmp1_masked, tmp2_masked);

    % Save results
    MAE_all(i)  = mae_val;
    RMSE_all(i) = rmse_val;
    SSIM_all(i) = ssim_val;
    PSNR_all(i) = psnr_val;

    % Single-frame results (optional)
    % disp(['Frame ', num2str(i), ...
    %       ': MAE=', num2str(mae_val), ...
    %       ', RMSE=', num2str(rmse_val), ...
    %       ', SSIM=', num2str(ssim_val), ...
    %       ', PSNR=', num2str(psnr_val)]);
end

% -------- Averages --------
MAE_mean  = mean(MAE_all,  'omitnan');
RMSE_mean = mean(RMSE_all, 'omitnan');
SSIM_mean = mean(SSIM_all, 'omitnan');
PSNR_mean = mean(PSNR_all, 'omitnan');

disp('=== Average metrics ===')
disp(['MAE mean  = ', num2str(MAE_mean)])
disp(['RMSE mean = ', num2str(RMSE_mean)])
disp(['SSIM mean = ', num2str(SSIM_mean)])
disp(['PSNR mean = ', num2str(PSNR_mean)])

%% Plot uncertainty distribution
i = 25;
test_data = abs(squeeze(matrix_sla_bayesnf(:,:,99,27))-squeeze(matrix_sla_bayesnf(:,:,99,1)));
test_data = test_data';
pro_3 = sla_dad(:,:,i);   
cloud_mask = isnan(pro_3);
test_data(~cloud_mask) = nan;  

% Visualize variance
figure;
m_pcolor(Lon_new,Lat_new,test_data);
m_gshhs('lc','patch',[0.5 0.5 0.5],'Edgecolor','k');
m_grid('xtick',[],'ytick',[],'tickdir','in','linest',':','xaxisloc','top');
colormap(flipud(nclCM(150,20)))
cb = colorbar;                 % Add colorbar
ylabel(cb,'\Delta SST/°C','Interpreter','tex');
cbarrow;
caxis([0,2])
title('Uncertainty Range (SST)'); % figure title
