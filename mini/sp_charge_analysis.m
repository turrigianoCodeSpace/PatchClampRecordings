%% Experiment names and folder pathes

%Name of the data folder
experiment = {'test'};

%location of the excised file
excised_data_fp = 'C:\Users\schum\Google_Drive\Lab\Data\culture_experiments\excised_data\';

%Filepath where you keep data folders
fp_data = 'C:\Users\schum\Google_Drive\Lab\Data\culture_experiments\mini\';

%whether to save results
save_results = 1;

%location to save the analyzed data
fp_analyzed_data = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\culture_experiments\analyzed_charge_results\';

%Name of the analyzed results
analyzed_file_name = strcat('CHARGE_ANALYSIS_',experiment{1},'.mat');

%sampling rate
sp_rate = 5000; %Hz

%plotting?
plot_on  = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HIPASS FILTER DESIGN %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Excised Data is run through a high pass filter to remove slow
%%%%%%% deviations from baseline

hfilt = designfilt('highpassiir', ...       % Response type
    'StopbandFrequency',0.35, ...     % Frequency constraints
    'PassbandFrequency',0.75, ...
    'PassbandRipple',0.5, ...
    'StopbandAttenuation',3, ...    % Magnitude constraints
    'SampleRate',5000);

%% Data readout (excised) and analysis
excised_data = load(strcat(excised_data_fp,'excised_',experiment{1,1},'.mat'));

excised_points = excised_data.excised_points;
aDAT = excised_data.data;
cell_id = excised_data.cell_id;
cell_num = numel(cell_id);

%pre-allocation
selected_Data = cell(1,cell_num);
TOTAL_TIME =cell(1,cell_num);
TRACE_CHARGE = cell(1,cell_num);


for ci = 1:cell_num
    if isempty(cell_id{1,ci})
        continue
    else
        trace_start = cell_id{1,ci}(1,1);
        trace_end = cell_id{1,ci}(end,1);

        for ti = trace_start:trace_end
            if ~ismember(ti,cell_id{1,ci})

               disp(strcat('Trace',num2str(ti), ' does not exist!'))
                TOTAL_TIME{1,ci}(:,ti) = NaN;
                TRACE_CHARGE{1,ci}(:,ti) = NaN;
                
                continue
            else
                pt_start = round(excised_points{1,1}{1,ci}{1,ti}(1,1));
                pt_end = round(excised_points{1,1}{1,ci}{1,ti}(2,1));
                
                %extract excised region and pass through a high-pass filter
                vals = aDAT{1,ci}(:,ti);
                vals_resampled(1:numel(resample(vals,1,2)),1) = resample(vals,1,2);
                if numel(vals_resampled) < pt_end
                    pt_end = numel(vals_resampled);
                end
                
                sDATA = filter(hfilt, detrend(vals_resampled(pt_start:pt_end)));
                sDATA = sDATA(isfinite(sDATA));
                selected_Data{1,ci}(1:numel(sDATA),ti) = sDATA;
                
                %calculate time and area under trace
                TOTAL_TIME{1,ci}(ti,1) = numel(sDATA)/sp_rate; %in seconds
                TRACE_CHARGE{1,ci}(ti,1) = abs(sum(sDATA)/sp_rate); %in pC
                
                if plot_on == 1
                    figure('position',[56 200 1400 490])
                    hold on
                    title(strcat('Cell',num2str(ci),': trace',num2str(ti)),'Interpreter', 'none')
                    plot(selected_Data{1,ci}(:,ti),'k')
                    
                    charge_str = strcat('charge = ', num2str(TRACE_CHARGE{1,ci}(ti,1)),' pC');
                    time_str = strcat('time = ', num2str(TOTAL_TIME{1,ci}(ti,1)),' s');
                    
                    str_dim = [0.65 0.6 0.3 0.3];
                    annotation('textbox',str_dim,'String',{charge_str, time_str},'FitBoxToText','on')
                end
                    
            end
        end
    end
end

%% save data

if save_results == 1
    
    cd(fp_analyzed_data)
    save(analyzed_file_name,'selected_Data','cell_id','TOTAL_TIME','TRACE_CHARGE')
end