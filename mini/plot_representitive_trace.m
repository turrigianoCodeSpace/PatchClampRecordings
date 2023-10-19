%% Select which cell and which trace to plot

% where to save the selected trace data
save_data_path = '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/representative_traces/mEPSC/';

% name of the saved data
save_data_name = 'APV_24h_1.mat';

% data folder
fp_data = '/Users/wwneuro/My_Drive/Lab/Data/culture_experiments/mini/';

experiment = '220902';

trace_file = 'cell8_0002.h5';

%% data readout and plotting
extracted_data_file = ws.loadDataFile(strcat(fp_data,experiment,'/',trace_file));
fields = fieldnames(extracted_data_file);  %get the fieldnames of the data struct
sprate = extracted_data_file.header.AcquisitionSampleRate;
selected_trace_data = extracted_data_file.(fields{2}).analogScans(:,1);

time_increment = 1/sprate;
number_of_datapoints = size(selected_trace_data,1); %get the number of data points of the sweep

%pick the interval you want to show (in seconds)
start_time = 45;
startpoint = start_time * sprate;
end_time = 45.6;
endpoint = end_time * sprate;
timepoints = (start_time:time_increment:end_time);

plot_data = detrend(selected_trace_data(startpoint:endpoint,1));

figure('position',[56 200 700 490]);
plot(timepoints,plot_data,'k-','Linewidth',1);
hold on
xlabel('Time (s)');
xlim([start_time end_time]);
ylim([-100 50])
    
%draw scale
plot([start_time+0.5; start_time+0.55],[-80; -80], '-k',[start_time+0.5;start_time+0.5],[-80; -70], '-k', 'LineWidth',2)
text(start_time+0.49, -76, '10 pA', 'HorizontalAlignment', 'right')
text(start_time+0.52, -84, '50 ms', 'HorizontalAlignment', 'center')
%set(gca, 'Visible', 'off')


%% save data
cd(save_data_path)
save(save_data_name, 'experiment','trace_file','start_time','end_time')