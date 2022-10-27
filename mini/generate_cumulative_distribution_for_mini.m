% This script loops through the all_mini_events_by_cell folder and groups
% them by experimental conditions: 
% For example:
    % DR+CNO- 1
    % CNO- 2
    % NT- 3

% The experimental condition of a cell can be uniquely determined by the date, 
% and the cell_num assigned on that day via a look-up cell_id_index table.  

% Once mini events have been grouped, 100 events will be randomly picked from
% each cell under each experimental condition, and then cumulative
% distributions will be generated for each condition, respectively.

%% Experimental settings
% in each ALL_EVENTS file, dates can be extracted from the last six digits

%save results
save_results = 1;

%experimental conditions (should match the numbers assigned to each group
%in the cell_id excel sheet)
exp_con = {'APV','Ctrl'};

%experiment name (correpsonding to the sheet name in the cell_id_index
% excel file)
exp_name = 'APV_24h_rise_1';

%where to save grouped files
fp_all_mini_group = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/all_mini_events_by_group/';

%location of ALL_EVENTS files
fp_all_events = '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/all_mini_events_by_cell/';

%subfolder (if any)
sub = 'rise_1';

%cell_id index readout
fp_mini = '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/mini_data_by_groups/';

cell_id_index = readtable(strcat(fp_mini,'cell_id_index.xlsx'),'Sheet',exp_name);


%% loop through the all events folder and categorize

%pre-allocation
all_mini_events = cell(1,numel(exp_con));

%loop through folder
cd(strcat(fp_all_events, sub))
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
    curr_date = curr_name(end-9:end-4);
    
    cy = strcat('20',curr_date(1:2));
    cM = curr_date(3:4);
    cday = curr_date(5:6);
    
    date = datetime(strcat(cy,'-',cM,'-',cday),'InputFormat','y-M-d',...
        'Format', 'M/d/y');
    if ~(ismember(date,cell_id_index.Date))
        continue
    else
        load(all_files(fi+2).name)
        
        for ci = 1:numel(all_amp)
            if isempty(all_amp{1,ci})
                continue
            else
                if isempty(find(cell_id_index.Date == date & cell_id_index.Cell_num == ci,1))
                    continue
                else
                
                    row = find(cell_id_index.Date == date & cell_id_index.Cell_num == ci);
                    cond_i = cell_id_index{row,'Cat'};
                    cell_ID = cell_id_index{row,'Cell_ID'};
                    all_mini_events{1,cond_i}{1,cell_ID} = all_amp{1,ci};
                end
            end
        end
    end
        
end

%% mini events selection
%pre-allocation
selected_mini_events = cell(1,numel(all_mini_events));

for ei = 1:numel(all_mini_events)
    
    counter = 1;
    
    %if less than 100 events in this cell, select all events
    %if more than 100 events, then randomly pick 100 events
    for ci = 1:numel(all_mini_events{1,ei})
        curr_cell = all_mini_events{1,ei}{1,ci};
        
        if size(curr_cell,1) < 100
            sel_events = curr_cell;
        else
            rand_index_exp = randi([1 size(curr_cell,1)],100,1);
            sel_events = curr_cell(rand_index_exp);
        end
        
        selected_mini_events{1,ei}(counter:counter+numel(sel_events)-1,1) = sel_events;
        counter = counter + numel(sel_events);
    end
end

%% generate cumulative statistics

%cumulative stats are saved in the cumulative cell array, in the order of
%the exp_con. For each condition, the first column stores amplitudes and
%the second one stores percentiles.

%pre-allocation
cumulative = cell(1,numel(selected_mini_events));

for eii = 1:numel(selected_mini_events)
    curr_cond = selected_mini_events{1,eii};
    
    min_amp = min(curr_cond);
    max_amp = max(curr_cond);
    
    [cu_X, cu_Y] = cumhist(curr_cond,[min_amp,max_amp],0.01);
    cumulative{1,eii}(:,1) = cu_X;
    cumulative{1,eii}(:,2) = cu_Y;
end

%% save to file

if save_results == 1
    cd(strcat(fp_all_mini_group, sub))
    save(strcat(exp_name,'.mat'),'all_mini_events','cell_id_index',...
        'selected_mini_events','cumulative','exp_con')
end