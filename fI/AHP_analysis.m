%%%% This script takes a analyzed f-I mat file and quantifies the
%%%% afterhyperpolarization following an action potential at the rheobase

%% change these for each run

% Name of file you would like to have data saved to
save_file_name = 'AHP_211215.mat';

% Name of the fI analyzed data file
fI_name = 'fI_211215.mat';

%filepath where you keep analyzed fI data
fp_fI = 'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\analyzed_fI_results\';

%location of the analyzed AHP file
fp_AHP = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\chronic_DREADDs\chronic_hm4di\analyzed_AHP';

%save results
save_results = 1;

%plot on
figure_on = 1;

%data selection mode: only applicable when plot_mode == 1
%1 for no normalization(start from the 1st current step)
%2 for normalization (start from the rheobase step)
select_mode = 2;

%stimulation amplitude: only applicable when plot_mode == 1
%indicates the stimulation amplitude (current injected, in pA)
stim = 0;

%% AHP quantification

% load analyzed fI file
cd(fp_fI)
load(fI_name)

% Pre-allocation
first_ap_data = cell(1,numel(cell_id));
first_ap_data_norm = cell(1,numel(cell_id));
first_ap_dV = cell(1,numel(cell_id));
AHP_peak_amp = NaN(numel(cell_id),1);
AHP_area = NaN(numel(cell_id),1);
AHP_duration = NaN(numel(cell_id,1));

for ci = 1:numel(cell_id)
    if isempty(aDAT{1,ci})
        continue
    else
        % Extract trace
        
        if select_mode == 1
            trace_id = cell_id{1,ci}(1,1) + stim/20 -1;
        elseif select_mode == 2
            trace_id = rheobase_ind(ci,1) + cell_id{1,ci}(1,1) + stim/20 - 1;
        end

        first_ap_start = V_th{1,ci}{1,trace_id}(1,1) - 60;
        
        if first_ap_start > 18000
            trace_id = trace_id + 1;
            first_ap_start = V_th{1,ci}{1,trace_id}(1,1) - 60;
        end

        
        data = aDAT{1,ci}(:,trace_id);

%         figure(1)
%         plot(data)

        % find the region of the first AP
        % starting from 6 ms before the threhold, use the first 5 ms for
        % baseline calculation

        if spike_count{1,ci}(trace_id,1) == 1
            first_ap_end = first_ap_start + 2000;
            
        else
            first_ap_end = V_th{1,ci}{1,trace_id}(2,1);
        end

        first_ap_vals = data(first_ap_start : first_ap_end);
        first_ap_dV_vals = diff(first_ap_vals);
        

%         figure(2)
%         plot(first_ap_vals)

        AHP_bl = nanmedian(first_ap_vals(1:50));
        first_ap_vals_norm = first_ap_vals - AHP_bl;

        [AHP_peak,peak_ind] = min(first_ap_vals_norm);
         

        % for high frequency trace, the duration might be an underestimate
        AHP_end_ind = find(first_ap_vals_norm(peak_ind:end) >= AHP_peak*0.1,1,'first') + peak_ind -1;
        if isempty(AHP_end_ind)
            [~,AHP_end_ind1] = max(first_ap_vals_norm(peak_ind:end));
            AHP_end_ind = AHP_end_ind1 + peak_ind -1;
        end

        AHP_duration_curr = (AHP_end_ind - peak_ind)/10; %in ms

        AHP_area_curr = sum(abs(first_ap_vals_norm(peak_ind : AHP_end_ind)));

        first_ap_data{1,ci} = first_ap_vals;
        first_ap_data_norm{1,ci} = first_ap_vals_norm;
        first_ap_dV{1,ci} = first_ap_dV_vals;
        AHP_peak_amp(ci,1) = AHP_peak;
        AHP_duration(ci,1)= AHP_duration_curr;
        AHP_area(ci,1) = AHP_area_curr;

        if figure_on == 1
            figure(ci)
            hold on
            plot(first_ap_vals_norm)
            for pi = peak_ind : AHP_end_ind
                plot([pi pi],[0 first_ap_vals_norm(pi)],'k-')
            end
        end
    end
end

%% save to file
if save_results == 1
    cd (fp_AHP)
    save(save_file_name,'first_ap_data','first_ap_data_norm','first_ap_dV','AHP_peak_amp','AHP_duration','AHP_area')
end