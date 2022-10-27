% This script takes grouped mini events, randomly picks 100 events from
% each cell under each experimental condition, and then calculates the
% cumulative distribution

%% load mini populations (grouped by experimental condition)

%whether to save results
save_results = 1;

%location of grouped mini events
fp_all_mini_group = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/all_mini_events_by_group/';

%sub folder
sub = 'rise_1';

%experiment name (correpsonding to the sheet name in the cell_id_index
% excel file)
exp_name = 'APV_24h_rise_1';

%load grouped mini events
cd(strcat(fp_all_mini_group, sub))
load(strcat(exp_name,'.mat'))

%where to save the analyzed cumulative data
fp_cumu = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/cumulative_analysis/';

%name of the saved cumulative data
save_file_name = strcat(exp_name, '_cumu_stats.mat');


%% prepare data
% choose experimental conditions you want to scale and plot for comparison
% refer to the exp_con for the order of cell arrays in cumulative
exp_ind = 1; %DR+CNO
ctrl_ind = 2; %CNO

%color code (darker analogous purple/blue)
color{1} = '#573E5C';
color{2} = '#913DA2'; 
color{3} = '#63B6FF';
color{4} = '#D095DB';


exp_amps = selected_mini_events{1,exp_ind}(:,1); 
ctrl_amps = selected_mini_events{1,ctrl_ind}(:,1); 

% make sure both groups have same number of selected mini events
% the group with more events is normalized to the one with less events
% (random selection)
if size(exp_amps,1) > size(ctrl_amps,1)
    rand_index = randi([1,size(exp_amps,1)],size(ctrl_amps,1),1);
    exp_amps_final = exp_amps(rand_index);
    ctrl_amps_final = ctrl_amps;
    
elseif size(exp_amps,1) < size(ctrl_amps,1)
    rand_index = randi([1,size(ctrl_amps,1)],size(exp_amps,1),1);
    exp_amps_final = exp_amps;
    ctrl_amps_final = ctrl_amps(rand_index);
    
else
    exp_amps_final = exp_amps;
    ctrl_amps_final = ctrl_amps;
end    

%sort event amplitudes from low to high
exp_amps_sort = sort(exp_amps_final);
ctrl_amps_sort = sort(ctrl_amps_final);

%% Scaling

%use linear fitting to get scale factor (scale exp to ctrl)
f = fittype('a*x+b','coefficients',{'a','b'});
slope_start = (max(exp_amps_sort)-min(exp_amps_sort))/(max(ctrl_amps_sort)-min(ctrl_amps_sort));
[lnfit,gof] = fit(ctrl_amps_sort,exp_amps_sort,f,'StartPoint',[slope_start 0]);
exp_scaled = (exp_amps_sort-lnfit.b) ./ lnfit.a;

%plot scaled data
figure()
plot(ctrl_amps_sort,ctrl_amps_sort,'k', 'Linewidth',2)
hold on
plot(ctrl_amps_sort,exp_amps_sort,'o','MarkerSize',12)
plot(lnfit)

ax1 = gca;
ax1.FontSize = 14;
ax1.LineWidth = 2;
ax1.YLabel.String = strcat(exp_con{exp_ind}, ' (pA)');
ax1.YLabel.FontSize = 14;
%ax1.XTick = [0 20 40 60 80];
ax1.XLabel.String = strcat(exp_con{ctrl_ind}, ' (pA)');
ax1.XLabel.FontSize = 14;
legend(strcat(exp_con{ctrl_ind}, ' vs. ', exp_con{ctrl_ind}),...
    strcat(exp_con{exp_ind}, ' vs. ', exp_con{ctrl_ind}),'FontSize',12,'Location','southeast')

legend('boxoff')
text(10,max(ctrl_amps_sort),strcat('y=',num2str(lnfit.a),'x',num2str(lnfit.b)),'FontSize',12)

box off

hold off

%% generate cumulative distribution plots

%cumulative distribution for scaled exp data
[cumu_exp_X_scaled, cumu_exp_Y_scaled] = cumhist(exp_scaled,[min(exp_scaled) max(exp_scaled)],0.01);

%two-sample Kolmogorov-Smirnov test
%after scaling
[h1,p1,ks2stat1] = kstest2(exp_scaled,ctrl_amps_sort);
%before scaling
[h2,p2,ks2stat2] = kstest2(exp_amps_sort,ctrl_amps_sort);


figure();
hold on

plot(cumulative{1,ctrl_ind}(:,1), cumulative{1,ctrl_ind}(:,2),'Color',color{ctrl_ind},'LineWidth',2)

%plot(cumu_exp_X_scaled, cumu_exp_Y_scaled,'k','LineWidth',2)

plot(cumulative{1,exp_ind}(:,1), cumulative{1,exp_ind}(:,2),'Color',color{exp_ind},'LineWidth',2)

ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;
ax.YLabel.String = 'Cumulative %';
ax.YLabel.FontSize = 14;
ax.XLim = [5 45];
ax.YLim = [0 110];
ax.YTick = [0 50 100];
ax.XLabel.String = 'mEPSC amplitude (pA)';
ax.XLabel.FontSize = 14;
legend(exp_con{ctrl_ind},'Scaled',exp_con{exp_ind},'FontSize',14,'Location','southeast')
legend('boxoff')
text(ax.XLim(2)-20,70,strcat(exp_con{exp_ind}, 'vs.',exp_con{ctrl_ind},': p=',num2str(p2)),'FontSize',12)
text(ax.XLim(2)-20,60,strcat('Scaled vs.',exp_con{ctrl_ind},': p=',num2str(p1)),'FontSize',12)
title(strcat(num2str(ax.XLim(2)),'pA cutoff'))
box off
hold off

%% save files
if save_results == 1
    cd(strcat(fp_cumu, rise))
    save(save_file_name,'exp_amps_sort','ctrl_amps_sort','exp_scaled','cumu_exp_X_scaled','cumu_exp_Y_scaled',...
        'h1','p1','ks2stat1','h2','p2','ks2stat2','rand_index')
end
