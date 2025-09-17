%% This program creates training and testing datasets based on CMEMS
%% 1. Create training and testing datasets for the original BayesNF
load("SST_dataset.mat");
start_date = datetime('2016-1-1'); % Start date
end_date = datetime('2016-2-28');  % End date
date_indices = find(time_sst >= start_date & time_sst <= end_date);
time = time_sst(date_indices);
sla_dad = temp_cloud(:,:,date_indices);
temp_gt = temp_GT(:,:,date_indices);
aim_zos = sla_dad;
aim_lon = Lon_new;
aim_lat = Lat_new;
output_1322 = [ones(size(aim_zos,3),1)];
excel_lon = aim_lon;
excel_lat = aim_lat;
num_1322 = hours((ti me-datetime(1900,1,1))/24+2/24)';

% for i = 1:size(aim_zos)
%     subplot(1,2,1)
%     pcolor(aim_lon,aim_lat,aim_zos(:,:,i))
%     shading flat
%     subplot(1,2,2)
%     pcolor(aim_lon,aim_lat,temp_gt(:,:,i))
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
excel_data(location_nan,:) = [];

%% 3. Create dataset
% train
train_dataset = excel_data;
num_rows = size(train_dataset, 1);   % Get the number of rows in the dataset
row_indices = (1:num_rows)';         % Create a sequence of row indices from 1 to num_rows
train_dataset = [row_indices, train_dataset];

% Set header
header = {'', 'location', 'datetime', 'longitude', 'latitude', 'sla'};
% Combine header and data
train_out = [header; num2cell(train_dataset)];
% Save as CSV
% excel_file = 'train_discussion_SST.csv';
excel_file = 'train_CMEMS_SST.csv';
writecell(train_out, excel_file);

%% validation
validation_dataset = excel_data_GT;
num_rows = size(validation_dataset, 1);
row_indices = (1:num_rows)';
validation_dataset = [row_indices, validation_dataset];

% Set header
header = {'', 'location', 'datetime', 'longitude', 'latitude', 'sla'};
% Combine header and data
validation_out = [header; num2cell(validation_dataset)];
% Save as CSV
% excel_file = 'validation_discussion_SST.csv';
excel_file = 'validation_CMEMS_SST.csv';
writecell(validation_out, excel_file);
