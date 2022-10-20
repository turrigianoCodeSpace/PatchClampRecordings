function [timepoints] = plot_trace(filename,ch,n)
%plot_trace('filename') plots the data as a function of time
%Inputs:   
%   filename- 'the complete name of the data file, including extensions'
%   (e.g. 'cell1_001.h5')
%   ch- channel ID (corresponds to the channel ID assignment on the
%   amplifier, ch = 1 under most cases if only one channel is active)
%   n- trace ID as shown in the filename. (e.g.
%   for'cell1_0003-0020.h5',n=3 returns trace 0003)
%   if n is omitted and the file includes multiple traces, then it
%   generates plots for all traces one by one upon pressing any key
%
%Outputs:
%   timepoints- x-axis converted into seconds (not necessary if only for plotting)

%Installation of DAQmax driver (from National Instruments) and WaveSurfer Required


if nargin<3 
    n=0;
end

extracted_data = ws.loadDataFile(filename);
time_increment = 1/extracted_data.header.AcquisitionSampleRate;
fields = fieldnames(extracted_data);

if extracted_data.header.NAOChannels == 1    
    y_unit = extracted_data.header.AIChannelUnits;
else
    y_unit = extracted_data.header.AIChannelUnits{ch,1};
end

x_unit = 'Time (s)';
figure_position = [56 200 1000 490];

total_filename = numel(filename);
trace_start = extracted_data.header.NextSweepIndex;

if n >0 
    n_norm = n-trace_start+1;
else
    n_norm = 0;
end


%convert x-axis unit to seconds
if n == 0 %no specific trace number indicated
    number_of_datapoints = size(extracted_data.(fields{n_norm+2}).analogScans(:,ch),1); 
else
    number_of_datapoints = size(extracted_data.(fields{n_norm+1}).analogScans(:,ch),1);
end

endtime = time_increment*number_of_datapoints;
timepoints = (0:time_increment:endtime);
timepoints = timepoints(1:end-1);

%locate specific trace or continuously display all traces in the file

if contains(filename,'-')
    if n_norm == 0
        for i = 1:(length(fields)-1)
            figure('position',figure_position);
            plot(timepoints,extracted_data.(fields{i+1}).analogScans(:,ch),'k-','Linewidth',1);
            xlabel(x_unit);
            xlim([0 endtime]);
            ylabel(y_unit);
            title(strcat(filename(1:total_filename-3),': ', fields(i+1)), 'Interpreter', 'none');

            waitforbuttonpress
        end
    else
    figure('position',figure_position);
    plot(timepoints,extracted_data.(fields{n_norm+1}).analogScans(:,ch),'k-','Linewidth',1);
    xlabel(x_unit);
    xlim([0 endtime]);
    ylabel(y_unit);
    title(strcat(filename(1:total_filename-3),': ', fields(n_norm+1)), 'Interpreter', 'none');
    
    end
else
    figure('position',figure_position);
    plot(timepoints,extracted_data.(fields{n_norm+2}).analogScans(:,ch),'k-','Linewidth',1);
    xlabel(x_unit);
    xlim([0 endtime]);
    ylabel(y_unit);
    title(strcat(filename(1:total_filename-3),': ', fields(n_norm+2)), 'Interpreter', 'none');
    
end

