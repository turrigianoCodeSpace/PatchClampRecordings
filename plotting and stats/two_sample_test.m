
%name of groups
groups = {'TTX_area','APV_area'};

%name of the variable 
parameter = 'MFR';

%Rheobase mode (values in pA or rheobase trace indicators)
rheo_mode = 2;

%testing mode: 0 for single column variable, 1 for multiple column
%variable
test_mode = 0;

%data selection mode: only applicable when plot_mode == 1
%1 for no normalization(start from the 1st current step)
%2 for normalization (start from the rheobase step)
select_mode = 1;

%stimulation amplitude: only applicable when plot_mode == 1
%indicates the stimulation amplitude (current injected, in pA)
stim = 300;

%stimulation increment (only applicable when plot_mode == 1)
%indicates the current injection increment during stimulation
inc = 25;

%confidence level
cl = 0.05;

sample = cell(1,numel(groups));
for gi = 1:numel(groups)
    if test_mode == 0
            sample{1,gi} = eval(strcat(groups{gi},'.',parameter));
        elseif test_mode ==1
            data_temp = eval(strcat(groups{gi},'.',parameter));
            if rheo_mode == 1
            rheo = eval(strcat(groups{gi},'.','Rheobase'));
            elseif rheo_mode == 2
            rheo = eval(strcat(groups{gi},'.','Rheobase_ind'));
            end
            sample{1,gi} = get_step_values(data_temp, inc, rheo, rheo_mode, select_mode, stim);
    end
end
%normality is tested using Shapiro-Wilk test
%if both samples follow normal distribution, both tests can be used
%otherwise the nonparametric MW test should be used
[h1,p1,W1] = swtest(sample{1,1});
[h2,p2,W2] = swtest(sample{1,2});

if h1 && h2 == 0
    %both are normal, two-sample t test
    [h,p,ci,stats] = ttest2(sample{1,1}, sample{1,2},'alpha',cl);
else
    %at least one sample does not follow normal distribution, Mann Whitney
    %U test
    [p,h,stats] = ranksum(sample{1,1}, sample{1,2},'alpha',cl);
end
