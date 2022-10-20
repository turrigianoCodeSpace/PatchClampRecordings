%This script records the X coordinates of raw mini data which will be used
%for analysis

%% to be changed for each new run

% specify experiment
experiment = {'test'};

% location to saved excised coordinates
fp_excised_pts = 'C:\Users\schum\Google_Drive\Lab\Data\culture_experiments\excised_data';

% location of mini data
fp_data = 'C:\Users\schum\Google_Drive\Lab\Data\culture_experiments\mini\';

%% start excising raw data

excised_points = cell(1,numel(experiment));
excised_filename = cell(1,numel(experiment));
abnormal_points = cell(1,numel(experiment));
abp = 0;

for jj = 1:numel(experiment)
    current_experiment = experiment{jj};
    current_folder = strcat(fp_data,experiment{1});
    excised_filename{1,jj} = strcat('excised_',experiment{1,jj},'.mat');
    
    [~, data, cell_id, cell_num, filename] = h5_file_readout(current_folder);
    
    file_num = numel(filename);
    excised_points{1,jj} = cell(1,file_num);
    
   for ci = 1:max(cell_num)
       if isempty(cell_id{1,ci}) == 1
           continue
       end
       disp(strcat('Current Cell: ', num2str(ci)))
       
       for ti = 1: size(cell_id{1,ci},1)
           
           trace_id = cell_id{1,ci}(ti,1);
           disp(strcat('Excising Trace', num2str(trace_id)))

            %datafile readout
            if isnan(data{1,ci}(1,trace_id))
                disp(strcat('Trace', num2str(trace_id), ' does not exist!'))
                continue
            else            
                vals = nonzeros(data{1,ci}(:,trace_id));               
            end

            %clear resampled
            clear resampled;
            %resample at 1/2
            resampled(1:numel(resample(vals,1,2)),1) = resample(vals,1,2); 
            %resample(x,p,q) resamples the input sequence x at p/q times
            %the original sample rate
            %resampled = vals; %remove the downsampling

            figure('position',[119 171 1210 611])
            hold on
            axis([1 numel(resampled) (mean(resampled,'omitnan')-120) (mean(resampled,'omitnan')+50)])
            %nanmean() returns the sample mean ignoring NaNs
            title(strcat('Cell',num2str(ci), ': trace', num2str(trace_id)),'Interpreter', 'none');
            plot(resampled)

            %disp(strcat(current_experiment,'_',num2str(g),' choose start and end points'))
            [X,~] = ginput; %ginput raises crosshairs in the current axes to identify points in the figure with the mouse
            excised_points{1,jj}{1,ci}{trace_id} = X;
            disp(strcat(num2str(numel(X)), ' points recorded for cell', num2str(ci), ' trace', num2str(trace_id)))
            
            if numel(X) ~= 2
                abp = abp + 1;
                abnormal_points{1,jj}(abp,1) = ci;
                abnormal_points{1,jj}(abp,2) = trace_id;
            end
            
            close all
            fclose all;
       end %trace
   end %cell
end %experiment

%% Reselect coordinates for traces having abnormal number of points (aka ~= 2)
reselect = input('Reselect? This will only modify traces that do not have 2 coordinates (1 = y, 2 = n):','s');
reselect = str2double(reselect);


abnormal_points_re = cell(1,numel(experiment));
abp = 0;

if reselect == 1
    if isempty(abnormal_points{1,jj})
        disp('All traces have been correctly truncated.')
    else
    
        for api = 1:size(abnormal_points{1,jj},1)
            ci_re = abnormal_points{1,jj}(api,1);
            ti_re = abnormal_points{1,jj}(api,2);
            vals_re = nonzeros(data{1,ci_re}(:,ti_re));  

            clear resampled;
            resampled(1:numel(resample(vals_re,1,2)),1) = resample(vals_re,1,2); 

            figure('position',[119 171 1210 611])
            hold on
            axis([1 numel(resampled) (mean(resampled,'omitnan')-120) (mean(resampled,'omitnan')+50)])
            title(strcat('Cell',num2str(ci_re), ': trace', num2str(ti_re)),'Interpreter', 'none');
            plot(resampled)
            
            [X,~] = ginput; 
            excised_points{1,jj}{1,ci_re}{ti_re} = X;
            disp(strcat(num2str(numel(X)), ' points recorded for cell', num2str(ci_re), ' trace', num2str(ti_re)))


            if numel(X) ~= 2
                abp = abp + 1;
                abnormal_points_re{1,jj}(abp,1) = ci_re;
                abnormal_points_re{1,jj}(abp,2) = ti_re;
            end
            
            close all
            fclose all;
        end
        
        abnormal_points = abnormal_points_re;
    end
end

%% Reselect coordinates for specfic traces
sp = input('Do you want to reselect coordinates for specifc traces?(1 = y, 2 = n):','s');

if str2num(sp) == 1
    
    cci = input('cell number:','s');
    cci = str2double(cci);

    tti = input('trace number:','s');
    tti = str2double(tti);

    vals_sp = nonzeros(data{1,cci}(:,tti));  

    clear resampled;
    resampled(1:numel(resample(vals_sp,1,2)),1) = resample(vals_sp,1,2); 

    figure('position',[119 171 1210 611])
    hold on
    axis([1 numel(resampled) (mean(resampled,'omitnan')-120) (mean(resampled,'omitnan')+50)])
    title(strcat('Cell',num2str(cci), ': trace', num2str(tti)),'Interpreter', 'none');
    plot(resampled)

    [X,~] = ginput; 
    excised_points{1,jj}{1,cci}{tti} = X;
    disp(strcat(num2str(numel(X)), ' points recorded for cell', num2str(cci), ' trace', num2str(tti)))

    close all
    fclose all;
end

%%
%%%%%%%%%%%%% SAVE  (careful with changes!)

result = input('Save results? This will overwrite previous file unless renamed! (1 = y, 2 = n):','s');

result = str2double(result);

if result == 1
    cd (fp_excised_pts)
       
    save(excised_filename{1},'data','excised_points','cell_id')

end
