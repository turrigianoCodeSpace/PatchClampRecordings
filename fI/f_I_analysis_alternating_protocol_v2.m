% This script is for data files that include both seal test and fI traces
% (i.e. in the h5 data file "cell1_0001-0020", all odd trace ids are seal 
% test, and all even traces are fI, or vice versa)
% All traces are read out from the h5 structure first, then saved
% separately in two cell structures: seal_test{} and f_I{}
% Passive and fI properties are then analyzed for each neuron,
% respectively.

%Should be noted that seal tests are also obtained under current clamp.

% Rheobase for each neuron is calculated from a data file that is different
% from the fI one; based on the first trace that shows AP in the fI file, 
% smaller current steps are used to determine the rheobase.
%% Experimental parameters setup (Always double check before each run)

%whether to save analyzed results (0 = no, 1 = yes)
save_results = 1; 

%filepath of fI data files
fp_fI_data = '/Users/wwneuro/My_Drive/Lab/Data/slice_NT/fI/';

%filepath of rheobase data files
fp_rheo_data = '/Users/wwneuro/My_Drive/Lab/Data/slice_NT/rheobase/';

%filepath of Vm data files
fp_vm_data = '/Users/wwneuro/My_Drive/Lab/Data/slice_NT/resting_vm/';

%filepath of analyzed fI data
fp_fI_analyzed_data = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/slice_NT/analyzed_fI_results';

%filepath of analyzed seal test data (aka passive properties)
fp_pp_analyzed_data = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/slice_NT/analyzed_seal_test';

%experiment(s) (defined as the 6-digit format of the recording date)
experiment = {'230831'};

%subfolders under a certain experiment if any (e.g. before/after)
sub = '';

%sampling rate (in Hz)
sample_rate = 10000;

%trace organization within the data file
%if the trace set starts with a seal test (i.e. seal tests are odd-numbered
%traces), then seal_odd = 1; otherwise seal_odd = 0;
seal_odd = 1; 

%%%% parameters for seal test pulse %%%%
%latency of pulse onset (in s)
step_start_pp = 0.5;

%duration of the seal test pulse (in seconds)
pulse_pp = 0.5;

%pulse amplitude (in A)
I_step = -50*10^-12;

%whether to show individual pp traces during analysis (1 = on, 0 = off)
figure_pp_on = 0;

%%%% parameters for fI/rheobase pulses
%latency of pulse onset (in s)
step_start_fI = 0.5;

%duration of the current injection pulse (in s)
pulse_fI = 1;

%current injection step amplitude (in pA)
curr_inc = 20;

%current injection step amplitude for rheobase (in pA)
curr_inc_rheo = 5;

%number of fI traces
%cti = 10;

%number of rheobase traces
%cti_rheo = 5;

%holding potential (in mV)
v_hold = -70;

%whether to show individual fI traces during analysis (1 = on, 0 = off)
figure_fI_on = 0;

%whether to show individual rheobase traces during analysis
figure_rheo_on = 0;

%% Data readout and categorization

for jj = 1:1%numel(experiment)
    current_experiment = experiment{jj};
    
    %file names of the analyzed fI and pp data were they to be saved
    save_name_fI = strcat('fI_', current_experiment, sub, '.mat'); %fI
    save_name_pp = strcat('pp_', current_experiment, sub, '.mat'); %pp
    
    current_fI_folder = strcat(fp_fI_data, current_experiment, '/', sub);
    [raw_h5_files, data, ~, cell_num, filename] = h5_file_readout(current_fI_folder);
    
                
    %pre-allocate cell arrays for fI and seal test traces
    aDAT_fI = cell(1, max(cell_num));
    aDAT_seal = cell(1, max(cell_num));
    %aDAT_vm = cell(1, max(cell_num));
    %aDAT_rheo = cell(1, max(cell_num));
    cell_id = cell(1, max(cell_num));
    %cell_id_rheo = cell(1, max(cell_num));
%%    
    %separate seal test and fI traces
    for ci = 1:max(cell_num)
        
        raw_trace_order = zeros(1,size(data{1,ci},2));
        %raw_trace_head = zeros(1,size(data{1,ci},2));
        %trace_break_ind = zeros(1,size(data{1,ci},2));
        
        for tii = 1:size(data{1,ci},2)
            if isnan(data{1,ci}(1,tii)) 
                aDAT_fI{1,ci}(1:size(data{1,ci},1),int8(tii/2)) = NaN;
                aDAT_seal{1,ci}(1:size(data{1,ci},1),int8(tii/2)) = NaN;
                cell_id{1,ci}(int8(tii/2),1) = 0; 
%                 if cell_id{1,ci}(int8(tii/2),1) > 0
%                     continue
%                 else
%                     cell_id{1,ci}(int8(tii/2),1) = 0; 
%                 end
            else
                
                raw_trace_order(1,tii) = tii;
                raw_trace_head = find(raw_trace_order,1,'first');
%                 if tii == 1
%                     continue
%                 else
%                     if tii > 1 && raw_trace_order(tii-1) == 0
%                         trace_break_ind(1,tii) = 1;
%                         raw_trace_head(1,tii) = raw_trace_head(find(raw_trace_head,1,'last'))+1;
% 
%                     else
%                         if sum(trace_break_ind)>0 && mod(sum(trace_break_ind),2) == 1
%                             raw_trace_head(1,tii) = raw_trace_head(1,tii-1);
%                         end
% 
%                     end
%                 end
                %raw_trace_order(1,tii) = tii;
                
                if mod(raw_trace_head,2) == 1 %trace ID in the raw file starts with an odd number
                    if mod(tii,2) == 0
                        if seal_odd == 1
                            aDAT_fI{1,ci}(1:numel(data{1,ci}(:,1)),tii/2) = ...,
                                data{1,ci}(1:numel(data{1,ci}(:,1)),tii);

                        else
                            aDAT_seal{1,ci}(1:numel(data{1,ci}(:,1)),tii/2) = ...,
                                data{1,ci}(1:numel(data{1,ci}(:,1)),tii);

                        end

                        cell_id{1,ci}(tii/2,1) = tii/2;
                    else
                        if seal_odd == 1
                            aDAT_seal{1,ci}(1:numel(data{1,ci}(:,1)),(tii-1)/2+1) = ...,
                                data{1,ci}(1:numel(data{1,ci}(:,1)),tii);

                        else
                            aDAT_fI{1,ci}(1:numel(data{1,ci}(:,1)),(tii-1)/2+1) = ...,
                                data{1,ci}(1:numel(data{1,ci}(:,1)),tii);

                        end

                        %cell_id{1,ci}((tii-1)/2+1,1) = (tii-1)/2+1;
                    end
                else  %trace ID in the raw file starts with an even number
                    if mod(tii,2) == 0 
                        if seal_odd == 1
                            aDAT_seal{1,ci}(1:numel(data{1,ci}(:,1)),tii/2) = ...,
                                data{1,ci}(1:numel(data{1,ci}(:,1)),tii);

                        else
                            aDAT_fI{1,ci}(1:numel(data{1,ci}(:,1)),tii/2) = ...,
                                data{1,ci}(1:numel(data{1,ci}(:,1)),tii);

                        end

                        cell_id{1,ci}(tii/2,1) = tii/2;
                    else
                        if seal_odd == 1
                            aDAT_fI{1,ci}(1:numel(data{1,ci}(:,1)),(tii-1)/2) = ...,
                                data{1,ci}(1:numel(data{1,ci}(:,1)),tii);

                        else
                            aDAT_seal{1,ci}(1:numel(data{1,ci}(:,1)),(tii-1)/2) = ...,
                                data{1,ci}(1:numel(data{1,ci}(:,1)),tii);

                        end

                        %cell_id{1,ci}((tii-1)/2,1) = (tii-1)/2;
                    end 
                end
            end
        end
        
        %Should have the same number of seal test and fI traces, if not,
        %discard the last trace.
        
        if mod(numel(nonzeros(raw_trace_order)),2) == 1
            cell_id{1,ci}(end,1) = 0;
            cell_id{1,ci} = nonzeros(cell_id{1,ci});
        end
        
    end
    
    %remove zeros from cell_id array
    for cci = 1:max(cell_num)
        cell_id{1,cci} = nonzeros(cell_id{1,cci});
    end
        
    %Rheobase data readout (if any)
    current_rheo_folder = strcat(fp_rheo_data, current_experiment, '/', sub);  
    
    if ~isfolder(current_rheo_folder)
        aDAT_rheo = NaN;
        cell_id_rheo = NaN;
    else
        [~, data_r, cell_id_rheo, cell_num_rheo, ~] = h5_file_readout(current_rheo_folder);
        aDAT_rheo = data_r;
    end
    
    %Vm data readout (if any)    
    current_vm_folder = strcat(fp_vm_data, current_experiment, '/', sub);
    
    if ~isfolder(current_vm_folder)
        aDAT_vm = NaN;
    else       
        [~, data_v, ~, ~, ~] = h5_file_readout(current_vm_folder);
        aDAT_vm = data_v;
    end
    
    %% Passive property calcualtion and firing characterization
    
    %%% pre-allocate data structures for pp analysis
    cell_stats_par_num = 7;
    Vm = cell(1,size(cell_id,2)); %resting membrane potential (no current inj)
    Vm_stats = NaN(max(cell_num),cell_stats_par_num);
    Rin = cell(1,size(cell_id,2)); %input resistance (in MOhm)
    Rin_stats = NaN(max(cell_num),cell_stats_par_num);
    Rin_outliers = cell(1,size(cell_id,2));
    Cm = cell(1,size(cell_id,2)); %membrane capacitance (in pF)
    Cm_stats = NaN(max(cell_num),cell_stats_par_num);
    Cm_outliers = cell(1,size(cell_id,2));
    fit_vals = cell(1,size(cell_id,2));
    fit_trace_time = cell(1,size(cell_id,2));
    exp_fit = cell(1,size(cell_id,2));
    cell_stats = cell(1,numel(experiment));
    burst_ind = NaN(max(cell_num),4);
    

    %%% pre-allocate data structures for fI analysis  
    MFR = cell(1,max(cell_num)); %mean firing rate
    ISI = cell(1,max(cell_num)); %interspike interval (for each two APs)
    IFR = cell(1,max(cell_num)); %instantenous firing rate (for each ISI)
    IFR_ave = cell(1,max(cell_num)); %average of IFRs for all spikes evoked in a single trace
    V_th = cell(1,max(cell_num)); %AP threshold (for each AP)
    V_th_ave = cell(1,max(cell_num)); %average AP threshold
    V_th_1st = cell(1,max(cell_num)); %threshold of the first AP in each trace
    V_th_4th = cell(1,max(cell_num)); %threshold of the fourth AP in each trace
    V_th_stats = NaN(max(cell_num),4); %stats of firing thresholds of the 1st APs 
    f_trough = cell(1,max(cell_num)); %AP fast trough (for each AP)
    f_trough_ave = cell(1,max(cell_num)); %average fast trough
    dV_sec = cell(1,max(cell_num)); %dV/dt for each trace
    AP_peak = cell(1,max(cell_num)); %peak amplitde of each AP
    AP_peak_ave = cell(1,max(cell_num)); %average AP peak amp
    AP_peak_1st = cell(1,max(cell_num)); %first AP amp
    Vrm = cell(1,max(cell_num)); %holding membrane potential
    Vrm_stat = NaN(max(cell_num),cell_stats_par_num); %stats of Vrm during current injection
    spike_count = cell(1,max(cell_num)); %number of spikes evoked
    curr_inj = cell(1,max(cell_num)); %current injected
    rheobase = NaN(max(cell_num),1); % each row stores the rheobase of a cell
    rheobase_ind = NaN(max(cell_num),1); %each row stores the index of rheobase (which stimulation)
    rheobase_sp_num = cell(1,max(cell_num)); %number of spikes detected in the rheobase data array for each cell
    lat = cell(1,max(cell_num)); %latency between the start of the current step and the first AP
    delta_isi = cell(1,max(cell_num)); %rate of change in interspike intevals
    adp_index = cell(1,max(cell_num)); %spike adaptation as the rate of ISI change of all ISIs within a trace
    adp_index_fl_ind = cell(1,max(cell_num)); %spike adaptation as the difference between the first and the last ISI
    udratio = cell(1,max(cell_num)); %upstroke/downstroke ratio
    udratio_slope = cell(1,max(cell_num)); %slope of the u/d ratio change across APs
    udratio_median = cell(1,max(cell_num)); %median of u/d ratio for each current step
    width = cell(1,max(cell_num)); %AP width at half height
    width_1st = cell(1,max(cell_num)); %width of the first AP for each current step
    width_slope = cell(1,max(cell_num)); %slope of the width change across APs
    width_median = cell(1,max(cell_num)); %median of width for each current step
    
    %% loop through all cells
    for ci =1:max(cell_num)
        if isempty(cell_id{1,ci}) == 1
            continue
        else
            trace_start = cell_id{1,ci}(1,1);
            trace_end = cell_id{1,ci}(end,1);

%             if size(cell_id{1,ci},1) > cti
%                 trace_end = cell_id{1,ci}(cti,1);
%             else
%                 trace_end = cell_id{1,ci}(end,1);
%             end
            
            for ti = trace_start:trace_end
                if ismember(ti,cell_id{1,ci}) == 0
                    Rin{1,ci}(ti,1) = NaN;
                    Cm{1,ci}(ti,1) = NaN;
                    fit_vals{1,ci}(ti,:) = NaN;
                    fit_trace_time{1,ci}(ti,:) = NaN;
                    exp_fit{1,ci}(ti,:) = NaN;
                    
                    MFR{1,ci}(ti,:) = NaN;
                    ISI{1,ci}(ti,:) = NaN;
                    IFR{1,ci}(ti,:) = NaN; 
                    V_th{1,ci}{ti} = {};
                    f_trough{1,ci}{ti} = {};
                    dV_sec{1,ci}{ti} = {};
                    AP_peak{1,ci}{ti} = NaN;
                    Vrm{1,ci}(ti,1) = NaN; 
                    spike_count{1,ci}(ti,1) = NaN;
                    curr_inj{1,ci}(ti,1) = NaN;
                    lat{1,ci}(ti,:) = NaN; 
                    adp_index{1,ci}(ti,:) = NaN;
                    udratio{1,ci}(ti,:) = NaN;
                    width{1,ci}(ti,:) = NaN;
                    delta_isi{1,ci}(ti,:) = NaN;

                    continue
                else

                    seal_data = aDAT_seal{1,ci}(:,ti);
                    [Rin_curr, Cm_curr, fit_vals_curr, trace_time_curr, exp_fit_curr] = ...,
                        get_PP_I_clamp(seal_data, step_start_pp, pulse_pp, I_step, sample_rate, figure_pp_on);
                    
                    Rin{1,ci}(ti,1) = Rin_curr;
                    Cm{1,ci}(ti,1) = Cm_curr;
                    fit_vals{1,ci}(1:numel(fit_vals_curr),ti) = fit_vals_curr;
                    fit_trace_time{1,ci}(1:numel(trace_time_curr),ti) = trace_time_curr;
                    exp_fit{1,ci}{ti} = exp_fit_curr;                    
                    
                    
                    fI_data = aDAT_fI{1,ci}(:,ti);
                    
                    %holding membrane potential (in mV)
                    Vrm{1,ci}(ti,1) = mean(fI_data(1000:2000),'omitnan');

                    %current injected (in pA)
                    
                    curr_inj{1,ci}(ti,1) = (ti-trace_start+1)*curr_inc;

                    %Detect AP in each trace
                    [dV_sec_temp, V_th_temp, f_trough_temp, sp_num_temp, udratio_temp] = ...,
                        get_vthresh(fI_data, step_start_fI, pulse_fI, sample_rate, figure_fI_on);


                    %number of spikes detected in each trace
                    spike_count{1,ci}(ti,1) = sp_num_temp;

                    % AP threshold
                    V_th{1,ci}{ti} = V_th_temp; %first col stores indices, second col stores potentials

                    %average AP threshold (start from the 4th AP)
                    V_th_ave{1,ci}(ti,1) = mean(V_th{1,ci}{1,ti}(4:end,2),'omitnan');
                    
                    %the first AP in each trace
                    V_th_1st{1,ci}(ti,1) = V_th{1,ci}{1,ti}(1,2);
                    
                    %the fourth AP in each trace
                    V_th_4th{1,ci}(ti,1) = V_th{1,ci}{1,ti}(4,2);


                    %AP fast trough
                    f_trough{1,ci}{ti} = f_trough_temp;%first col stores indices, second col stores potentials

                    %average fast trough (all APs)
                    f_trough_ave{1,ci}(ti,1) = mean(f_trough{1,ci}{1,ti}(:,2),'omitnan');

                    %dV/dt for each trace
                    dV_sec{1,ci}{ti} = dV_sec_temp;

                    %MFR (mean firing rate)
                    MFR{1,ci}(ti,1) = sp_num_temp/pulse_fI; %mean FR in Hz (duration 1s)

                    %ISI (interspike interval, in ms)
                    %defined as the difference between the fast trough of the
                    %initial AP and the threshold of the next AP
                    if sp_num_temp == 0 || sp_num_temp == 1
                        ISI{1,ci}(ti,:) = 0;
                    else
                        for ii = 1:(sp_num_temp-1)
                            ISI{1,ci}(ti,ii) = (V_th_temp(ii+1,1) - f_trough_temp(ii,1))*0.1;
                        end
                    end

                    %IFR (instantaneous firing rate, in Hz))
                    %calculated as the reciprocal of ISI
                    if sp_num_temp == 0 || sp_num_temp == 1
                        IFR{1,ci}(ti,:) = 0;
                    else
                        for fi = 1:(sp_num_temp-1)
                            IFR{1,ci}(ti,fi) = 1000/ISI{1,ci}(ti,fi);
                        end

                    end
                    
                    %AP peak amplitude
                    %first col stores amplitudes, second col stores peak indice
                    if sp_num_temp == 0
                        AP_peak{1,ci}{ti} = NaN;
                    else
                        for pi = 1:sp_num_temp
                            [peak_amp,peak_ind] = ...
                                max(fI_data(int32(V_th_temp(pi,1)):int32(f_trough_temp(pi,1))));
                            AP_peak{1,ci}{1,ti}(pi,1) = peak_ind +V_th_temp(pi,1); %peak index
                            AP_peak{1,ci}{1,ti}(pi,2) = peak_amp; %peak amplitude
                        end
                    end

                    %calulate average/1st AP peak amplitude for each trace
                    if ~isnan(AP_peak{1,ci}{1,ti})
                        AP_peak_ave{1,ci}(ti,1) = mean(AP_peak{1,ci}{1,ti}(:,2),'omitnan');
                        AP_peak_1st{1,ci}(ti,1) = AP_peak{1,ci}{1,ti}(1,2);
                    else
                        AP_peak_ave{1,ci}(ti,1) = NaN;
                        AP_peak_1st{1,ci}(ti,1) = NaN;
                    end

                    %latency (in ms)
                    if sp_num_temp == 0
                       lat{1,ci}(ti,1) = NaN;
                    else
                       lat{1,ci}(ti,1) = (V_th_temp(1,1)-step_start_fI*sample_rate)*0.1;
                    end

                    %spike adaptation index
                    %method 1: calculated as the normalized difference
                    %between the last and the first ISI
                    if sp_num_temp < 3
                        adp_index_fl_ind{1,ci}(ti,1) = 0;
                    else
                        ISI_temp = nonzeros(ISI{1,ci}(ti,:));
                        adp_index_fl_ind{1,ci}(ti,1) = (ISI_temp(end)-ISI_temp(1))/ISI_temp(1); %should > 1
                    end
                    
                    
                    %the rate of ISI change
                    if sp_num_temp < 3
                        %adp_index{1,ci}(ti,1) = NaN;
                        delta_isi{1,ci}(ti,:) = NaN;
                        adp_index{1,ci}(ti,1) = 0;
                    else
                        %adp_index{1,ci}(ti,1) = nanmean(diff(ISI{1,ci}(ti,:)))/nansum(ISI{1,ci}(ti,:));
                        for adpi = 2:sp_num_temp-1
                            delta_isi{1,ci}(ti,adpi) = (ISI{1,ci}(ti,adpi) - ISI{1,ci}(ti,adpi-1))/...,
                                (ISI{1,ci}(ti,adpi) + ISI{1,ci}(ti,adpi-1));
                        end
                        
                        adp_index{1,ci}(ti,1) = sum(delta_isi{1,ci}(ti,:),'omitnan')/(sp_num_temp-2);
                    end
                    
%                     %spike adaptation index
%                     %the rate of ISI change
%                     if sp_num_temp < 3
%                         adp_index{1,ci}(ti,1) = NaN;
%                     else
%                         adp_index{1,ci}(ti,1) = mean(diff(ISI{1,ci}(ti,:)),'omitnan')...,
%                             /sum(ISI{1,ci}(ti,:),'omitnan');
%                     end

                    %upstroke/downstroke ratio
                    %defined as the ratio of the maxmium to the minimum dV/dt
                    %for each AP
                    if sp_num_temp == 0
                        udratio{1,ci}(ti,1:sp_num_temp) = NaN;
                        udratio_median{1,ci}(ti,1) = NaN;
                    else
                        udratio{1,ci}(ti,1:sp_num_temp) = udratio_temp(1:sp_num_temp);

                        %calculate median
                        udratio_median{1,ci}(ti,1) = median(udratio{1,ci}(ti,1:sp_num_temp), 'omitnan');

                        %calculate the slope of change in U/D ratio for each trace
                        if sp_num_temp < 2
                            udratio_slope{1,ci}(ti,1) = NaN;
                        else
                            [fm1,gof] = fit((1:sp_num_temp)',(udratio{1,ci}(ti,1:sp_num_temp))','poly1');

    %                         subplot(4,5,ti-trace_start+1);
    %                         plot((1:sp_num_temp)',(udratio{1,ci}(ti,1:sp_num_temp))','o')
    %                         hold on
    %                         plot(fm1)

                            udratio_slope{1,ci}(ti,1) = fm1.p1;
                        end
                    end


                    %AP width at half height (in ms)
                    if sp_num_temp == 0
                        width{1,ci}(ti,:) = NaN;
                        width_1st{1,ci}(ti,1) = NaN;
                        width_median{1,ci}(ti,1) = NaN;
                    else
                        for wi = 1:sp_num_temp
                            %downward halfheight
                            halfheight_d = 0.5*(AP_peak{1,ci}{1,ti}(wi,2)-f_trough{1,ci}{1,ti}(wi,2))+...
                                f_trough{1,ci}{1,ti}(wi,2);
                            halfheight_d_ind = find(fI_data(int32(AP_peak{1,ci}{1,ti}(wi,1)):int32(f_trough{1,ci}{1,ti}(wi,1)))<=...
                                halfheight_d,1,'first')+AP_peak{1,ci}{1,ti}(wi,1);

                            %extrapolate to the upward side
                            data_extrp = fI_data(int32(V_th{1,ci}{1,ti}(wi,1)):int32(AP_peak{1,ci}{1,ti}(wi,1)));
                            [data_extrp, index] = unique(data_extrp);
                            %data_extrp = data_extrp';
                            ind_extrp = V_th{1,ci}{1,ti}(wi,1):AP_peak{1,ci}{1,ti}(wi,1);
                            ind_extrp = ind_extrp(index)';


                            if isnan(interp1(data_extrp, ind_extrp, halfheight_d))
                               warning(strcat('Cannot find width at half-height for one or more AP(s) in Cell', ...
                                   num2str(ci), ' Trace', num2str(ti), '!'))
%                                 data_extrp = fI_data((int32(V_th{1,ci}{1,ti}(wi,1))-50)...,
%                                     :int32(AP_peak{1,ci}{1,ti}(wi,1)));
%                                 warning(strcat('Cell ', ci, 'Trace ', ti, ' has APs with high threshold!'))
%                                 [data_extrp, index] = unique(data_extrp);
%                                 ind_extrp = (V_th{1,ci}{1,ti}(wi,1)- 50):AP_peak{1,ci}{1,ti}(wi,1);
%                                 ind_extrp = ind_extrp(index)';
%                                 halfheight_u_ind = interp1(data_extrp, ind_extrp, halfheight_d);
                            else
                                halfheight_u_ind = interp1(data_extrp, ind_extrp, halfheight_d);
                            end
    %                         %upward halfheight
    %                         halfheight_u = 0.5*AP_peak{1,ci}(ti,wi)-V_th{1,ci}{1,ti}(wi,1);
    %                         half_height_u_ind = find(data(V_th{1,ci}{1,ti}(wi,2):AP_peak{1,ci}{1,ti}(wi,2))>=...
    %                             AP_peak{1,ci}{1,ti}(wi,1)-halfheight_u,1,'first');

                            %calculate width
                            width{1,ci}(ti,wi) = (halfheight_d_ind - halfheight_u_ind)*0.1;

                            %troubleshooting width

%                             data_ind = ind_extrp(1):int32(f_trough{1,ci}{1,ti}(wi,1));
%                             data = fI_data(data_ind);
%                             figure;
%                             hold on;
%                             plot(data_ind,data,'k')
%                             plot(halfheight_d_ind,halfheight_d,'o','markersize',10,'markerfacecolor','r')
%                             plot(halfheight_u_ind,halfheight_d,'o','markersize',10,'markerfacecolor','b')


                        end

                        %width of the first AP in each current step
                        width_1st{1,ci}(ti,1) = width{1,ci}(ti,1);
                        
                        %calculate median of width
                        width_median{1,ci}(ti,1) = median(width{1,ci}(ti,:),'omitnan');

                        %calculate width slope   
                        if sp_num_temp < 2
                            width_slope{1,ci}(ti,1) = NaN;
                        else
                            [fm2,gof] = fit((1:sp_num_temp)',(width{1,ci}(ti,1:sp_num_temp))','poly1');

    %                         subplot(4,5,ti-trace_start+1);
    %                         plot((1:sp_num_temp)',(width{1,ci}(ti,1:sp_num_temp))','o')
    %                         hold on
    %                         plot(fm2)

                            width_slope{1,ci}(ti,1) = fm2.p1;
                        end
                    end
                    
                    %mean IFR (mean of all IFRs)
                    if sp_num_temp > 1
                        IFR_ave{1,ci}(ti,1) = mean(nonzeros(IFR{1,ci}(ti,:)),'omitnan');
                    else
                        IFR_ave{1,ci}(ti,1) = 0;
                    end                  
                end %non-empty trace
            end %per trace           
        end 
        
        %Burst activity identification (probably only applies to cortical
        %bursting PYNs that only have 2-3 spikes in a burst)
        %Definition: at the first spike train that contains 3-4 APs (3 ISIs),
        %if 2 APs in a burst: 2nd ISI > 3* 1st ISI
        %if 3 APs in a burst: 3rd ISI > 3* average of the first two ISIs
        %mean ISI during the burst <= 30 ms
        
                              
        if isempty(find(spike_count{1,ci} >= 3,1,'first'))
            burst_ind(ci,1) = 0; % burst indicator, 1 = Y 0 = N
            burst_ind(ci,2) = NaN; % # of spikes in a burst
            burst_ind(ci,3) = NaN; % average ISI during the burst (in ms)
            burst_ind(ci,4) = NaN;
        else
             bi = find(spike_count{1,ci} >= 3,1,'first');

            if ISI{1,ci}(bi,1)*3 < ISI{1,ci}(bi,2) && ISI{1,ci}(bi,1) <= 30
                burst_ind(ci,1) = 1; % burst indicator, 1 = Y 0 = N
                burst_ind(ci,2) = 2; % # of spikes in a burst
                burst_ind(ci,3) = ISI{1,ci}(bi,1); % average ISI during the burst (in ms)
                burst_ind(ci,4) = ISI{1,ci}(bi,2);
            else
                if bi+1 > cell_id{1,ci}(end,1) || spike_count{1,ci}(bi+1,1)<4
                    burst_ind(ci,1) = 0; % burst indicator, 1 = Y 0 = N
                    burst_ind(ci,2) = NaN; % # of spikes in a burst
                    burst_ind(ci,3) = ISI{1,ci}(bi,1); % average ISI during the burst (in ms)
                    burst_ind(ci,4) = ISI{1,ci}(bi,2);
                else
                    if ISI{1,ci}(bi+1,1)*3 < ISI{1,ci}(bi+1,2) && ISI{1,ci}(bi+1,1) <= 30
                        burst_ind(ci,1) = 1; % burst indicator, 1 = Y 0 = N
                        burst_ind(ci,2) = 2; % # of spikes in a burst
                        burst_ind(ci,3) = ISI{1,ci}(bi+1,1); % average ISI during the burst (in ms)
                        burst_ind(ci,4) = ISI{1,ci}(bi+1,2);
                    else

                        if mean(ISI{1,ci}(bi,1:2),'omitnan')*3 < ISI{1,ci}(bi,3) && mean(ISI{1,ci}(bi,1:2),'omitnan') <= 30
                            burst_ind(ci,1) = 1; % burst indicator, 1 = Y 0 = N
                            burst_ind(ci,2) = 3; % # of spikes in a burst
                            burst_ind(ci,3) = mean(ISI{1,ci}(bi,1:2),'omitnan'); % average ISI during the burst (in ms)
                            burst_ind(ci,4) = ISI{1,ci}(bi,3);
                        else
                            burst_ind(ci,1) = 0; % burst indicator, 1 = Y 0 = N
                            burst_ind(ci,2) = NaN; 
                            burst_ind(ci,3) = mean(ISI{1,ci}(bi,1:2),'omitnan');
                            burst_ind(ci,4) = ISI{1,ci}(bi,3);
                        end
                    end
                end
            end
        end

        
        %calculate rheobase for each cell
        
        for tti = 1:(trace_start+size(cell_id{1,ci},1)-1)

            %if spike_count{1,ci}(tti,1) ~= 0 && spike_count{1,ci}(tti+1,1) ~= 0
            if spike_count{1,ci}(tti,1) ~= 0
                if tti == trace_start+size(cell_id{1,ci},1)-1               
                    rheo_est = curr_inj{1,ci}(tti,1);
                    rheobase_ind(ci,1) = tti-trace_start+1;    
                else
                    if spike_count{1,ci}(tti+1,1) ~= 0
                        rheo_est = curr_inj{1,ci}(tti,1);
                        rheobase_ind(ci,1) = tti-trace_start+1;    
                        break
                    end
                end
     
            else
                rheo_est = NaN;
                rheobase_ind(ci,1) = NaN;
            end
        end
        
        %get the exact rheobase from the rheobase data array
        if ~iscell(cell_id_rheo) || ~ismember(ci, cell_num_rheo)
            rheobase(ci,1) = rheo_est;
        else
            for rheo_i = 1:size(cell_id_rheo{1,ci},1)
                if ci > size(aDAT_rheo,2) || isempty(cell_id_rheo{1,ci}) 
                    continue
                else
                    rheo_data = aDAT_rheo{1,ci}(:,cell_id_rheo{1,ci}(rheo_i,1));

                    [~,~,~,sp_num_rheo,~] = get_vthresh(rheo_data,...
                                         step_start_fI, pulse_fI,sample_rate,figure_rheo_on);
                    rheobase_sp_num{1,ci}(rheo_i,1) = sp_num_rheo;

                end
            end


            %if no AP detected in the rheobase data array, use the rheobase
            %estimated from the fI data array
            if ci > size(aDAT_rheo,2) || isempty(cell_id_rheo{1,ci}) || ~sum(rheobase_sp_num{1,ci}(:,1), 'omitnan')
                rheobase(ci,1) = rheo_est;
            else
                delta_rheo = find(rheobase_sp_num{1,ci}(:,1)', 1, 'first')*curr_inc_rheo;

                if (rheo_est-curr_inc+delta_rheo) > rheo_est
                    rheobase(ci,1) = rheo_est;
                else               
                    rheobase(ci,1) = rheo_est-curr_inc + delta_rheo;
                end

            end
        end
        
        %calculate Vm from Vm data array (several 30s traces were taken before
        %the fI interrogations)
        if ~iscell(aDAT_vm)
            Vm{1,ci} = NaN;
            Vm_stats(ci,1:4) = NaN;          
        else
            
             if ci > size(aDAT_vm,2) || isempty(aDAT_vm{1,ci})
                 Vm{1,ci} = NaN;
                 Vm_stats(ci,1:4) = NaN;  
             else
                 
                for vmi = 1:size(aDAT_vm{1,ci},2)
                    vm_data_temp = aDAT_vm{1,ci}(:,vmi);
                    vm_data_temp = vm_data_temp(5*sample_rate:25*sample_rate); %use the region between 5s and 25s
                    vm_mean = mean(vm_data_temp,'omitnan');
                    vm_std = std(vm_data_temp,'omitnan');

                    for vmii = 1:size(vm_data_temp,1)
                        if (vm_data_temp(vmii,1)-vm_mean) > 3*vm_std
                            vm_data_temp(vmii,1) = NaN;
                        end
                    end

                    Vm{1,ci}(vmi,1) = mean(vm_data_temp,'omitnan');
                end

                if size(Vm{1,ci},1)<3
                    if abs(diff(Vm{1,ci}(:,1)))/abs(Vm{1,ci}(end,1)) > 0.1
                        Vm_stats(ci,1) = Vm{1,ci}(end,1);
                    else
                        Vm_stats(ci,1) = mean(Vm{1,ci}(:,1),'omitnan');
                    end
                else

                    if std(Vm{1,ci},'omitnan') / mean(Vm{1,ci},'omitnan') > 0.1
                        Vm_stats(ci,1) = Vm{1,ci}(end,1);
                    else
                        Vm_stats(ci,1) = mean(Vm{1,ci}(:,1),'omitnan');
                    end
                end

                if Vm_stats(ci,1) > -50 % cells with Vm > -50mV should be excluded (in 2nd col.)
                    Vm_stats(ci,2) = 0;
                else
                    Vm_stats(ci,2) = 1;
                end
                
             end   
        end
        
        
        %mean and std of Vrm during current injections for each cell
        Vrm_stat(ci,1) = mean(Vrm{1,ci}(trace_start:ti,1),'omitnan');%1st col- average
        Vrm_stat(ci,2) = std(Vrm{1,ci}(trace_start:ti,1),'omitnan');%2nd col- std
        Vrm_stat(ci,3) = -(Vrm_stat(ci,2)/Vrm_stat(ci,1))*100; %coefficient of variation
        Vrm_stat(ci,4) = abs(Vrm_stat(ci,1)-v_hold)/abs(v_hold)*100; %deviation from V_hold 
        %whether Vrm is acceptable (CV should be less than 10%, deviation
        %from V_hold should be less than 5%)
        if Vrm_stat(ci,4)>5 || Vrm_stat(ci,3) > 10
            Vrm_stat(ci,5) = 0;
        else
            Vrm_stat(ci,5) = 1;
        end
        
        %Threshold test
        %Cells whose average firing threshold (1st AP of each trace) exceeds -25mV should be excluded
        V_th_stats(ci,1) = mean(V_th_1st{1,ci}(trace_start:trace_end),'omitnan');%1st col- average
        V_th_stats(ci,2) = std(V_th_1st{1,ci}(trace_start:trace_end),'omitnan');%2nd col- std
        V_th_stats(ci,3) = -(V_th_stats(ci,2)/V_th_stats(ci,1))*100; %coefficient of variation
        
        if V_th_stats(ci,1) > -30
            V_th_stats(ci,4) = 0;
        else
            V_th_stats(ci,4) = 1;
        end
                
        %Rin and Cm stats (CV should be less than 15% after outlier removal (if any))
        Rin{1,ci} = nonzeros(Rin{1,ci});
        Rin_stats(ci,1) = mean(Rin{1,ci}(:,1),'omitnan'); %1st col- average
        Rin_stats(ci,2) = std(Rin{1,ci}(:,1),'omitnan'); %2nd col- std
        Rin_stats(ci,3) = (Rin_stats(ci,2)/Rin_stats(ci,1))*100; %coefficient of variation
        
        Cm{1,ci} = nonzeros(Cm{1,ci});
        Cm_stats(ci,1) = mean(Cm{1,ci}(:,1),'omitnan');%1st col- average
        Cm_stats(ci,2) = std(Cm{1,ci}(:,1),'omitnan');%2nd col- std
        Cm_stats(ci,3) = (Cm_stats(ci,2)/Cm_stats(ci,1))*100; %coefficient of variation
        
        if Rin_stats(ci,3) > 15
            if sum(isoutlier(Rin{1,ci}(:,1),'median')) > 0 %outlier detected by median absoulte deviations
                out_ind = isoutlier(Rin{1,ci}(:,1),'median');
                Rin_rev = Rin{1,ci}(~out_ind,1);
                Rin_outliers{1,ci} = out_ind;              
                Rin_stats(ci,4) = sum(out_ind);
                
                if Rin_stats(ci,4) < 3
                    Rin_stats(ci,5) = mean(Rin_rev,'omitnan'); %mean after outlier removal
                    Rin_stats(ci,6) = std(Rin_rev,'omitnan'); %std after outlier removal
                    Rin_stats(ci,7) = Rin_stats(ci,6)/Rin_stats(ci,5)*100; %revised CV
                    
                    
                    if Rin_stats(ci,7) < 15
                        Rin_stats(ci,4) = 1;
                    else
                        Rin_stats(ci,4) = 0;
                    end
                else
                    Rin_stats(ci,4) = 0;
                    
                end
            else   
                Rin_stats(ci,4) = 0;
                out_ind = isoutlier(Rin{1,ci}(:,1),'median');
                Rin_outliers{1,ci} = out_ind;
            end
                               
        else
            Rin_stats(ci,4) = 1;
            
        end
        
        if Cm_stats(ci,3) > 15
            if sum(isoutlier(Cm{1,ci}(:,1),'median')) > 0 %outlier detected by median absoulte deviations
                out_ind = isoutlier(Cm{1,ci}(:,1),'median');
                Cm_rev = Cm{1,ci}(~out_ind,1);
                Cm_outliers{1,ci} = out_ind;              
                Cm_stats(ci,4) = sum(out_ind);
                
                %sum(out_ind)
                if Cm_stats(ci,4) < 3
                    Cm_stats(ci,5) = mean(Cm_rev,'omitnan'); %mean after outlier removal
                    Cm_stats(ci,6) = std(Cm_rev,'omitnan'); %std after outlier removal
                    Cm_stats(ci,7) = Cm_stats(ci,6)/Cm_stats(ci,5)*100; %revised CV
                
                    if Cm_stats(ci,7) < 15
                        Cm_stats(ci,4) = 1;
                    else
                        Cm_stats(ci,4) = 0;
                    end
                else
                    Cm_stats(ci,4) = 0; %exclude if more than two outliers detected
                end
            else
                Cm_stats(ci,4) = 0; % if no outlier detected but CV too high
            end
                               
        else
            Cm_stats(ci,4) = 1;
        end
        
        
    end %per cell
    
        %All cell stats (Rin, Cm, Vrm, etc) should be saved in a separate
        %data structure
        %outliers in a separate data structure
        outliers.Rin = Rin_outliers;
        outliers.Cm = Cm_outliers;
        
        cell_stats{1,jj}.Rin_stats = Rin_stats;
        cell_stats{1,jj}.Cm_stats = Cm_stats;
        cell_stats{1,jj}.Vrm_stats = Vrm_stat;
        cell_stats{1,jj}.V_th_stats = V_th_stats;
        cell_stats{1,jj}.Vm_stats = Vm_stats;
        cell_stats{1,jj}.outliers = outliers;
        cell_stats{1,jj}.burst_ind = burst_ind;

end %per experiment

%% Save results
if save_results

    %%% save a copy for passive properties only
    cd(fp_pp_analyzed_data)
    save(save_name_pp, 'aDAT_seal','aDAT_vm','Vm','Vm_stats','Rin','Cm','Rin_stats',...,
        'Cm_stats','fit_vals','fit_trace_time','exp_fit')

    %%% save everything (PP and FI) in the analyzed fI file
    cd(fp_fI_analyzed_data)
    save(save_name_fI, 'aDAT_fI','aDAT_rheo','cell_id','cell_id_rheo','cell_stats','adp_index',...
        'adp_index_fl_ind','AP_peak','AP_peak_1st','AP_peak_ave','Cm','curr_inj','delta_isi','raw_h5_files',...
        'dV_sec','experiment','f_trough','f_trough_ave','IFR','IFR_ave','ISI','lat','MFR','rheobase',...
        'rheobase_ind','rheobase_sp_num','Rin','spike_count','udratio','udratio_median','udratio_slope',...
        'v_hold','V_th','V_th_1st','V_th_4th','V_th_ave','Vrm','Vm','width','width_1st','width_median','width_slope');

end