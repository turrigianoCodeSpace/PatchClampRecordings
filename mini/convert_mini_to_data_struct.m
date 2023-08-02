% This script loops through the analyzed_mini_results + waveform_average_kinetics folder and groups
% them by experimental conditions: 
% For example: (these indice should match 'Cat' in the cell_id excel file) 
    % Ctrl- 1
    % APV- 2
    % TTX- 3

% The experimental condition of a cell can be uniquely determined by the date, 
% and the cell_num assigned on that day via a look-up cell_id_index table.  

%% %% Change these accordingly based on how you want to group the data
% in each anaylzed_mini_results file, dates can be extracted from the last six digits

%save results
save_results = 1;

%rise time cutoff
rise_cutoff = 'rise_1';

%experiment name (correpsonding to the sheet name in the cell_id_index
% excel file)
exp_name = strcat('TTX_GLYX13_24h_',rise_cutoff);

%file name
saved_file_name = 'TTX_GLYX13_24h.mat';

%where to save grouped files
fp_grouped_data = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/mini_data_by_groups/';

%experimental conditions
exp_con = {'Ctrl','TTX','TTX_GLYX'};

%import cell_id_index table 
cd('/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/mini_data_by_groups')
cell_id_index = readtable('cell_id_index.xlsx','Sheet',exp_name);

%location of analyzed mini results
fp_analyzed_mini = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/analyzed_mini_results/';

%location of waveform average kinetics results
fp_wavg_kinetics = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/waveform_average_kinetics/';

%% extract and group parameters from analyzed mini files
data_struct_temp = cell(1,numel(exp_con));

cd(strcat(fp_analyzed_mini,rise_cutoff))
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
    if ismember('m',curr_name(14:19))
        continue
    else
        curr_date = curr_name(14:19);
    end
    
    
    cy = strcat('20',curr_date(1:2));
    cM = curr_date(3:4);
    cday = curr_date(5:6);
    
    date = datetime(strcat(cy,'-',cM,'-',cday),'InputFormat','y-M-d',...
        'Format', 'M/d/y');
    
    if ~(ismember(date,cell_id_index.Date))
        continue
    else
        load(all_files(fi+2).name)
        
        for ci = 1:size(cell_id,2)
            if isempty(cell_id{1,ci})
                continue
            else
                if isempty(find(cell_id_index.Date == date & cell_id_index.Cell_num == ci,1))
                    continue
                else
                
                    row = find(cell_id_index.Date == date & cell_id_index.Cell_num == ci);
                    cond_i = cell_id_index{row,'Cat'};
                    cell_ID = cell_id_index{row,'Cell_ID'};
                    ti_start = cell_id{1,ci}(1,1);
                    trace_ct = size(cell_id{1,ci},1);
                    ti_end = cell_id{1,ci}(trace_ct,1);
                    
                    data_struct_temp{cond_i}.amp(cell_ID,1) = cell_stats{1,1}.amp_stats(ci,1);
                    data_struct_temp{cond_i}.frq(cell_ID,1) = cell_stats{1,1}.frq(ci,1);
                    
                    %for passive properties
                    if isnan(cell_stats{1,1}.Rin_stats(ci,5))
                        data_struct_temp{cond_i}.Rin(cell_ID,1) = cell_stats{1,1}.Rin_stats(ci,1);
                    else
                        data_struct_temp{cond_i}.Rin(cell_ID,1) = cell_stats{1,1}.Rin_stats(ci,5);
                    end
                    
                    if isnan(cell_stats{1,1}.Cm_stats(ci,5))
                        data_struct_temp{cond_i}.Cm(cell_ID,1) = cell_stats{1,1}.Cm_stats(ci,1);
                    else
                        data_struct_temp{cond_i}.Cm(cell_ID,1) = cell_stats{1,1}.Cm_stats(ci,5);
                    end
                    
                    if isnan(cell_stats{1,1}.Vm_stats(ci,5))
                        data_struct_temp{cond_i}.Vm(cell_ID,1) = cell_stats{1,1}.Vm_stats(ci,1);
                    else
                        data_struct_temp{cond_i}.Vm(cell_ID,1) = cell_stats{1,1}.Vm_stats(ci,5);
                    end
                    
                    if isnan(cell_stats{1,1}.Rs_stats(ci,5))
                        data_struct_temp{cond_i}.Rs(cell_ID,1) = cell_stats{1,1}.Rs_stats(ci,1);
                    else
                        data_struct_temp{cond_i}.Rs(cell_ID,1) = cell_stats{1,1}.Rs_stats(ci,5);
                    end
                    
                end
            end
        end
    end
end




%% extract and group parameters from waveform average kinetics files

clear all_files
clear cell_id

cd(strcat(fp_wavg_kinetics,rise_cutoff))
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
    if ismember('m',curr_name(1:6))
        continue
    else
        curr_date = curr_name(1:6);
    end
    
    
    cy = strcat('20',curr_date(1:2));
    cM = curr_date(3:4);
    cday = curr_date(5:6);
    
    date = datetime(strcat(cy,'-',cM,'-',cday),'InputFormat','y-M-d',...
        'Format', 'M/d/y');
    
    if ~(ismember(date,cell_id_index.Date))
        continue
    else
        load(all_files(fi+2).name)
        
        for ci = 1:size(cell_id,2)
            if isempty(cell_id{1,ci})
                continue
            else
                if isempty(find(cell_id_index.Date == date & cell_id_index.Cell_num == ci,1))
                    continue
                else
                
                    row = find(cell_id_index.Date == date & cell_id_index.Cell_num == ci);
                    cond_i = cell_id_index{row,'Cat'};
                    cell_ID = cell_id_index{row,'Cell_ID'};
                    ti_start = cell_id{1,ci}(1,1);
                    trace_ct = size(cell_id{1,ci},1);
                    ti_end = cell_id{1,ci}(trace_ct,1);
                    
                    data_struct_temp{cond_i}.risetime(cell_ID,1) = riseAve{1,1}(ci,1);
                    data_struct_temp{cond_i}.decaytau_1st(cell_ID,1) = decayTau_1{1,1}(ci,1);
                    data_struct_temp{cond_i}.decaytau_2nd(cell_ID,1:2) = decayTau_2{1,1}(ci,1:2);
                    data_struct_temp{cond_i}.decaytot_1st(cell_ID,1) = decayTot_1{1,1}(ci,1);
                    data_struct_temp{cond_i}.decaytot_2nd(cell_ID,1:2) = decayTot_2{1,1}(ci,1:2);
                    data_struct_temp{cond_i}.decay_2nd_comp_prct(cell_ID,1:2) = component_percentage{1,1}(ci,1:2);
                    data_struct_temp{cond_i}.charge(cell_ID,1) = wavg_charge{1,1}(ci,1);
                                       
                end
            end
        end
    end
end


%% CHANGE FOR EACH EXPERIMENT BEFORE SAVING
Ctrl = data_struct_temp{1,1};
TTX = data_struct_temp{1,2};
TTX_GLYX = data_struct_temp{1,3};

if save_results == 1
    
    cd(fp_grouped_data)
    save(saved_file_name,exp_con{1,:},'exp_con','cell_id_index')
end