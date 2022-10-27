%%%% This script extracts a specified action potential from all
%%%% cells for each experiment (e.g. first AP at rheobase)

%% Experimental settings

%save results
save_results = 1;

%where to save the extracted AP data
fp_fap = '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/firstAP';

%location of the analyzed f-I files
fp_fI = '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/analyzed_fI_results/';

%name of the experiment
exp = '220923';

%data selection mode:
%1 for no normalization(start from the 1st current step)
%2 for normalization (start from the rheobase step)
select_mode = 2;

%Rheobase mode (values in pA or rheobase trace indicators)
rheo_mode = 2;

%stimulation amplitude: only applicable when plot_mode == 1
%indicates the stimulation amplitude (current injected, in pA)
stim = 0;

%stimulation increment
%indicates the current injection increment during stimulation
inc = 25;

%%
% load fI file
file_name = strcat('firstAP_', exp, '.mat');
fI_name = strcat('fI_',exp,'.mat');
load(strcat(fp_fI,fI_name))


stim_ind = NaN(numel(cell_id),1);
trace_ind = NaN(numel(cell_id),1);
V_firstAP = NaN(150,numel(cell_id));
dV_firstAP = NaN(150,numel(cell_id));

V_firstAP_ave = NaN(150,1);
dV_firstAP_ave = NaN(150,1);

trace_length = NaN(30,1);
width = NaN(30,1);

for si = 1:numel(cell_id)
   if isnan(rheobase_ind(si,1))
       continue
   else
        
       if select_mode == 1
           if stim == 0
               warning('Trace 0 does not exist!')
           else
               stim_ind(si,1) = stim/inc;
           end
       elseif select_mode == 2
           stim_ind(si,1) = stim/inc + rheobase_ind(si,1); 
       else
           warning('Unknown mode')
       end

        trace_ind(si,1) = stim_ind(si,1) + cell_id{1,si}(1,1) - 1;
        trace_current = trace_ind(si,1);
        
        V_th_temp = V_th{1,si}{1,trace_current}(1,1);
        f_trough_temp = f_trough{1,si}{1,trace_current}(1,1);
        trace_start = V_th_temp - 10;
        trace_end = f_trough_temp + 10;
        
        %trace length: threshold-10:f_trough+10
        trace_length(si,1) = trace_end - trace_start;

        V_firstap = aDAT_fI{1,si}((trace_start : trace_end),trace_current);

        V_firstAP(1:numel(V_firstap),si) = V_firstap;

        dV_firstap = dV_sec{1,si}{1,trace_current}((V_th{1,si}{1,trace_current}(1,1)-10)...,
            :(f_trough{1,si}{1,trace_current}(1,1)+10),1);
        dV_firstAP(1:numel(V_firstap),si) = dV_firstap;
   end
%     figure
%     plot(V_firstap)
%     hold on
%     plot(dV_firstap./10)
%     hold off
%     
end



% figure
% hold on
% for sii = 1:numel(cell_id)
%     plot(V_firstAP(:,sii))
%     plot(V_firstAP(:,sii), dV_firstAP(:,sii))
% end

% hold off

% figure
% plot(V_firstAP_ave)
% hold on
% plot(dV_firstAP_ave./10)
% hold off

% figure
% plot(V_firstAP_ave,dV_firstAP_ave)

%% save data
if save_results == 1
    cd(fp_fap)
    save(file_name,'V_firstAP','dV_firstAP', 'trace_length')
end