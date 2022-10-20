%%%% This script takes the imported f-I mat file and calculate basic stats
%%%% for each property

%%%% Under each field, data from each experimental condition (i.e. WT,
%%%% WT_CNO, DR_CNO, and, WT_saline), are saved in corresponding columns 1, 2, 3,
%%%% and 4.

%% Initialization

%location where the mat file will be saved
fp_analyzed_data = '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/fI_data_by_groups';

%name of the mat file
filename = {'ttx_apv_co_16h_stats.mat'};

%experimental conditions (stats in each column follow this order)
cond = {'Ctrl','APV','TTX'};

%functions for calculating stats
func = {@mean, @std, @nansem};

%field names (corresponding to properties)
fields = {'MFR','IFR','mean_IFR'};
%{'MFR', 'IFR', 'mean_IFR', 'Threshold', 'ADP_ind', 'Rheobase',...
%    'Lat', 'Rin', 'Cm', 'Udratio', 'Width_med','Width_first'};

%plot
plot_on = 1;

%save results
save_results = 0;

%% calculate stats 

%calculate stats for each current step
struct_all = cell(1,3); %results from each function are saved in corresponding cells

for sti = 1:numel(func) % per stat
    for cdi = 1:numel(cond) % per experimental condition
        curr_c = eval(cond{cdi});
        
        for fi = 1:numel(fields) % per field
            if size(curr_c.(fields{fi}),2)<2 % ignore fields that are not current-step related
                continue
            else
                for cui = 1:size(curr_c.(fields{fi}),1) % per current step
                    struct_all{sti}.(fields{fi})(cui,cdi) = func{sti}(curr_c.(fields{fi})(cui,:));
                end
            end
        end
    end
end

%current injection vector (just pick one vector, should all be the same in
%terms of plotting
current_inj = curr_c.curr_inj(:,1);

%transfer results to corresponding data structs for storage
ave = struct_all{1};
std = struct_all{2};
sem = struct_all{3};

%calculate median for each cell (med)
% and
%get values from the rheobase step for each cell (rheo)
%get values from 100 pA over rheobase for each cell (above100)

for cdii = 1:numel(cond) % per experimental condition
        curr_c = eval(cond{cdii});
        
        for fii = 1:numel(fields) % per field
            if size(curr_c.(fields{fii}),2)<2 % ignore fields that are not current-step related
                continue
            else
                for cuii = 1:size(curr_c.(fields{fii}),2) % per cell
                    med.(fields{fii})(cuii,cdii) = median(curr_c.(fields{fii})(:,cuii),'omitnan');
                end
                
%                 vals = get_step_values(curr_c.(fields{fii}), 2, 0);
%                 rheo.(fields{fii})(1:size(vals,1),cdii) = vals;
%                 
%                 vals_1 = get_step_values(curr_c.(fields{fii}), 2, 100);
%                 above100.(fields{fii})(1:size(vals_1,1),cdii) = vals_1;
                
            end
        end
end
 
%current injection vector (just pick one vector, should all be the same in
%terms of plotting

%% plotting

% data range for plotting
plot_range = 1:10;
plot_variable = 'IFR';
plot_stat = 'ave';
err = 'sem';
y_label = 'IFR (Hz)';
mksz = 12; %marker size

%tint factors
tint_factor = [0 0 0 0];
%tint and convert function
%x-RGB value, y-tint factor
fun = @(x,y) (x+(255-x)*y)./255;

%color code
% color2 = fun([255 203 164],tint_factor(2)); %deep peach
% color3 = fun([153 186 221],tint_factor(3)); %carolina blue
%color3 = fun([255 150 0],tint_factor(3)); %dark orange
% color1 = fun([250 128 114],tint_factor(1)); %salmon
%color4 = fun([70 130 180],tint_factor(4)); %steel blue
% 
%color2 = fun([255 203 164],tint_factor(3)); %deep peach
%color4 = fun([153 186 221],tint_factor(4)); %carolina blue

%CT color scheme
% color1 = '#DB7093'; %pink
% color2 = '#4680B2'; %steel blue

%CPP color scheme
% color1 = '#F17F73'; %salmon
% color2 = '#4680B2'; %steel blue
% color3 = '#F9A36B'; %light orange
% color4 = '#63B6FF'; %light blue

%analogous purple
% color1 = '#E984F9'; 
% color2 = '#849EF9';
% color3 = '#AD8DF2'; 

%darker color transition (still purple)
color1 = '#D095DB'; 
color2 = '#573E5C';
color3 = '#913DA2'; 
%color4 = '#F17F73'; %salmon
%color4 = '#63B6FF'; %light blue


curr_y = eval(strcat(plot_stat,'.',plot_variable));
curr_y_err = eval(strcat(err,'.',plot_variable));

if plot_on == 1

    figure
    hold on
    %cond1
    errorbar(current_inj(plot_range),curr_y(plot_range,1),curr_y_err(plot_range,1),'-o','MarkerSize',mksz,...
            'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'LineWidth',2,'Color',color1)
    %cond2
    errorbar(current_inj(plot_range),curr_y(plot_range,2),curr_y_err(plot_range,2),'-o','MarkerSize',mksz,...
            'MarkerEdgeColor',color2, 'MarkerFaceColor',color2,'LineWidth',2,'Color',color2)    
        
    %cond3
    errorbar(current_inj(plot_range),curr_y(plot_range,3),curr_y_err(plot_range,3),':^','MarkerSize',mksz,...
            'MarkerEdgeColor',color3, 'MarkerFaceColor','none','LineWidth',2,'Color',color3)
%     %cond4
%     errorbar(current_inj(plot_range),curr_y(plot_range,4),curr_y_err(plot_range,4),':^','MarkerSize',mksz,...
%             'MarkerEdgeColor',color4, 'MarkerFaceColor','none','LineWidth',2,'Color',color4)
        
%     %cond5
%     errorbar(current_inj(plot_range),curr_y(plot_range,5),curr_y_err(plot_range,5),':^','MarkerSize',mksz,...
%             'MarkerEdgeColor',color5, 'MarkerFaceColor','none','LineWidth',2,'Color',color5)
    %WT+saline
%     errorbar(current_inj(plot_range),curr_y(plot_range,2),curr_y_err(plot_range,2),'-o','MarkerSize',8,...
%             'MarkerEdgeColor','k', 'MarkerFaceColor','k','LineWidth',2,...
%             'Color','k')

    ax = gca;
    ax.FontSize = 14;
    ax.LineWidth = 2;
    ax.YLabel.String = y_label;
    ax.YLabel.FontSize = 14;
    ax.XLabel.String = 'Injected Current (pA)';
    ax.YLabel.FontSize = 14;
    ax.YLabel.Interpreter = 'none';
    ax.XLim = [0 260];
    ax.YLim = [0 Inf];

    legend('Ctrl','APV','TTX','FontSize',14,'Location','northwest','Interpreter','none','Box','off')
    %title('Bursting')

end

%% save results
if save_results == 1
    cd (fp_analyzed_data)
    save(filename{1},'ave','std','sem','med')
end
