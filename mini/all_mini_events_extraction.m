% This script extract all amplitudes of mini events from an analyzed mini data struct and 
% save events from each cell separately.

%% 
%whether to save results
save_results = 1;

%location of mini analysis files
fp_analyzed_mini = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/analyzed_mini_results/';

%subfolders (if any)
sub = 'rise_1';

%experiment
exp = '220903';

%location to save the ALL_EVENTS file
fp_all_events = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/all_mini_events_by_cell/';


%% Amplitude extraction

%name of the the saved file
mini_file_name = strcat('MINIANALYSIS_', exp, '.mat');
save_file_name = strcat('ALL_EVENTS_', exp,'.mat');

%load file
curr_exp = load(strcat(fp_analyzed_mini, sub, '/', mini_file_name));
curr_exp_amp = curr_exp.AMP_ALL;

cell_num = numel(curr_exp_amp{1,1});

%pre-allocation
all_amp_raw = cell(1,cell_num);
all_amp = cell(1,cell_num);

%extract amplitudes from each trace
tracect = 1;

for ci = 1:cell_num
    for ti = 1:size(curr_exp_amp{1,1}{1,ci},2)
        if sum(~isnan(curr_exp_amp{1,1}{1,ci}(:,ti))) == 0 || ...
                mean(curr_exp_amp{1,1}{1,ci}(:,ti),'omitnan') == 0
            continue
        else       
            lgth = size(curr_exp_amp{1,1}{1,ci}(:,ti),1);
            all_amp_raw{1,ci}(1:lgth,tracect) = curr_exp_amp{1,1}{1,ci}(:,ti);
            tracect = tracect + 1;
        end
    end
end

%convert NaNs into zeros, and then remove zeros before concatnenating all
%columns into a single column

for cii = 1:cell_num
    nan_ind = ~isnan(all_amp_raw{1,cii});
    all_amp{1,cii} = nonzeros(all_amp_raw{1,cii}(nan_ind));

end


%% save results
if save_results == 1
    cd(strcat(fp_all_events, sub))
    save(save_file_name,'all_amp')
end