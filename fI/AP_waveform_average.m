% This script loops through the firstAP folder and groups
% them by experimental conditions: 
% For example:
    % Ctrl- 1
    % APV- 2
    % TTX- 3

% The experimental condition of a cell can be uniquely determined by the date, 
% and the cell_num assigned on that day via a look-up cell_id_index table.  

% Once the AP waveform of all cells have been grouped by their
% experimental conditions, three types of mega waveform average will be
% generated: from unscaled aps, aps scaled to the median peak amp within
% the condition, and aps scaled to the max peak amp within the condition.
% Choose which waveform average to plot at the end of the script.

%% Experimental settings

%save results
save_results = 1;

%figure on
figure_on = 1;

%experiment name (correpsonding to the sheet name in the cell_id_index
% excel file)
exp_name = 'TTX_APV_final_dataset';

%where to save grouped files
fp_ap_group = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/fI_data_by_groups/firstAP_by_group';

%experimental conditions
exp_con = {'TTX','APV','TTX_APV','Ctrl_combined'};

%import cell_id_index table 
cd('/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/fI_data_by_groups')
cell_id_index = readtable('cell_id_index.xlsx','Sheet',exp_name);

%location of single action potential files (per cell)
AP_fp = '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/firstAP/';

%% group data by condition

%pre-allocation
all_ap = cell(1,numel(exp_con));
all_ap_norm = cell(1,numel(exp_con));

%loop through folder
cd(AP_fp)
all_files_temp = dir;
for fi = 1:size(all_files_temp,1)
    if strcmp('.DS_Store',all_files_temp(fi).name)
        delete '.DS_Store'
        continue
    end
end

all_files = dir;
file_num = numel(all_files)-2;

for fi = 1:file_num
    curr_name = all_files(fi+2).name;
    curr_date = curr_name(9:14);
    
    cy = strcat('20',curr_date(1:2));
    cM = curr_date(3:4);
    cday = curr_date(5:6);
    
    date = datetime(strcat(cy,'-',cM,'-',cday),'InputFormat','y-M-d',...
        'Format', 'M/d/y');
    if ~(ismember(date,cell_id_index.Date))
        continue
    else
        load(all_files(fi+2).name)
        
        for ci = 1:size(V_firstAP,2)
            if sum(V_firstAP(:,ci)) == 0
                continue
            else
                if isempty(find(cell_id_index.Date == date & cell_id_index.Cell_num == ci,1))
                    continue
                else
                
                    row = find(cell_id_index.Date == date & cell_id_index.Cell_num == ci);
                    cond_i = cell_id_index{row,'Cat'};
                    cell_ID = cell_id_index{row,'Cell_ID'};
                    all_ap{1,cond_i}(:,cell_ID) = V_firstAP(:,ci);
                    
                    bl = mean(all_ap{1,cond_i}(1:5,cell_ID));
                    all_ap_norm{1,cond_i}(:,cell_ID) = all_ap{1,cond_i}(:,cell_ID) - bl;
                end
            end
        end
    end
        
end


%% calculate meta average (all normalized to baseline)
%pre-allocation
all_meta_ap = cell(1,numel(exp_con)); %average
all_meta_ap_med = cell(1,numel(exp_con));%peak_scaled_average (to the median)
all_meta_ap_max = cell(1,numel(exp_con));%peak_scaled_average (to the maximum)

all_ap_med = cell(1,numel(exp_con)); %peak scaled to the median
all_ap_max = cell(1,numel(exp_con)); %peak scaled to the maximum
scale_med = cell(1,numel(exp_con));
scale_max = cell(1,numel(exp_con));

%calculate average for each row
for ei = 1:numel(exp_con)
    max_amp = NaN(size(all_ap{1,ei},2),1);
    scale = NaN(size(all_ap{1,ei},2),1);

        %peak-scaled
    for ci = 1:size(all_ap{1,ei},2)
        [max_val,~] = max(all_ap_norm{1,ei}(:,ci));
        max_amp(ci,1) = max_val;
    end

    for cii = 1:size(all_ap{1,ei},2)
        %scale peak amp of all cells to the median
        amp_med = median(max_amp,'omitnan');
        scale_med{1,ei}(cii,1) = amp_med/max_amp(cii,1);
        all_ap_med{1,ei}(:,cii) = all_ap_norm{1,ei}(:,cii) * scale_med{1,ei}(cii,1);

        [amp_max,~] = max(max_amp);
        scale_max{1,ei}(cii,1) = amp_max/max_amp(cii,1);
        all_ap_max{1,ei}(:,cii) = all_ap_norm{1,ei}(:,cii) * scale_max{1,ei}(cii,1);
    end

    for ri = 1:size(all_ap{1,ei},1)
        %not peak-scaled
        all_meta_ap{1,ei}(ri,1) = mean(all_ap_norm{1,ei}(ri,:),'omitnan');

        %peak-scaled
        all_meta_ap_med{1,ei}(ri,1) = mean(all_ap_med{1,ei}(ri,:),'omitnan');
        all_meta_ap_max{1,ei}(ri,1) = mean(all_ap_max{1,ei}(ri,:),'omitnan');
    end

end

%% plot

%Choose which AP average to plot
ap_plot = all_meta_ap_max;

%color code (darker analogous purple/blue)
color{4} = '#D095DB'; 
color{1} = '#573E5C';
color{2} = '#913DA2'; 
color{3} = '#63B6FF'; %light blue

plot_range = (1:110);

if figure_on == 1
    figure();
    hold on

    for gi = 1:numel(exp_con)

        plot(ap_plot{1,gi}(plot_range),'Color',color{gi},'LineWidth',4)

    end

    title('First AP at rheobase step')
    legend(exp_con{1,:},'Interpreter','none','Box','off')
    
    %draw scale
    plot([60; 80],[30; 30], '-k',[60;60],[30; 40], '-k', 'LineWidth',2)
    text(60, 35, '10 mV', 'HorizontalAlignment', 'right')
    text(70, 25, '20 ms', 'HorizontalAlignment', 'center')
    
    box off
    hold off
end

%% save to file

if save_results == 1
    cd(fp_ap_group)
    save(strcat(exp_name,'.mat'),'all_ap','all_ap_norm','all_ap_med','all_ap_max',...,
        'all_meta_ap','all_meta_ap_med','all_meta_ap_max','cell_id_index','exp_con')
end