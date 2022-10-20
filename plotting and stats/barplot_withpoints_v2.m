%% groups for plotting
%name of groups
groups = {'Ctrl_area','TTX_area','APV_area'};

%label for each group
g_label = {'Ctrl','TTX','APV'};

%name of the variable 
parameter = 'MFR';

%name of the varibale for label
parlabel = 'MFR area';

%unit of the variable
unit = 'pA*Hz';

%scale for unit conversion
scale = 1;

%exponent value of the Y axis
y_axis_expo = 0;

%lower and upper limit of the Y axis
%y_limit = [-inf inf]; %inf indicates auto [-inf,inf]

%direction of the Y axis
y_dir_all = {'normal','reverse'};
y_dir = y_dir_all{1};

%Y axis lower limit (0 for origin, Inf for automatically calcuated limit)
y_lo_lim = Inf;

%plotting mode: 0 for single column variable, 1 for multiple column
%variable
plot_mode = 0;

%data selection mode: only applicable when plot_mode == 1
%1 for no normalization(start from the 1st current step)
%2 for normalization (start from the rheobase step)
select_mode = 1;

%Rheobase mode (values in pA or rheobase trace indicators)
rheo_mode = 2;

%stimulation amplitude: only applicable when plot_mode == 1
%indicates the stimulation amplitude (current injected, in pA)
stim = 150;

%stimulation increment (only applicable when plot_mode == 1)
%indicates the current injection increment during stimulation
inc = 25;

%offset of scatter points with respect to the error bar
dp_offset = 0.2;

%whether to plot data in boxplots
box_on = 0;

%whether to compare medians (only if box_on == 1)
% on- show notched boxplots
% off- show regular boxplots
notch_on = 'off';

%tint factors (btw 0 and 1, closer to 1 = more tint)
tint_factor_1 = [0 0 0 0 0]; %bar tint
tint_factor_2 = [0 0 0 0 0]; %point tint

%tint and convert function
    %x-RGB value, y-tint factor
fun = @(x,y) (x+(255-x)*y)./255;

%col1-3 store face color, col4-6 store edge color
%bar

col_pp1 = 1:3;
col_pp2 = 4:6;

% col_bar = NaN(numel(groups),3);
% col_bar(1,col_pp1) = fun([153 186 221],tint_factor_1(1)); %carolina blue
% col_bar(2,col_pp1) = fun([70 130 180],tint_factor_1(2)); %steel blue
% col_bar(3,col_pp1) = fun([255 150 0],tint_factor_1(3)); %dark orange
% %col_bar(4,col_pp1) = fun([250 128 114],tint_factor_1(3)); %salmon
% 
% col_bar(1,col_pp2) = fun([153 186 221],tint_factor_2(1)); %carolina blue
% col_bar(2,col_pp2) = fun([70 130 180],tint_factor_2(2)); %steel blue
% col_bar(3,col_pp2) = fun([255 150 0],tint_factor_2(3)); %dark orange
% %col_bar(4,col_pp2) = fun([250 128 114],tint_factor_1(3)); %salmon


col_bar = cell(1,numel(groups));

% %deep peach and carolina blue
% col_bar{2} = '#FFCBA4';
%col_bar{3} = '#99BADD';

%salmon and steel blue
% col_bar{1} = '#FA8072';
%col_bar{2} = '#4682B4';
%col_bar{3} = '#FF9600'; %dark orange

%pink and blue (ctail color scheme)
% col_bar{1} = '#DB7093';
% col_bar{2} = '#4682B5';

% CPP color scheme
% col_bar{1} = '#F17F73'; %salmon
% col_bar{2} = '#4680B2'; %steel blue
% col_bar{3} = '#F9A36B'; %light orange
% col_bar{4} = '#63B6FF'; %light blue

%culture data color scheme
%darker color transition (still purple)
col_bar{1} = '#D095DB'; 
col_bar{3} = '#913DA2';
col_bar{2} = '#573E5C';
%col_bar{4} = '#F17F73'; %salmon
%col_bar{4} = '#63B6FF'; %light blue

%points
col_poi = NaN(numel(groups),3);
%col_poi(1,col_pp1) = fun([0 0 0],tint_factor_2(1)); %black
% col_poi(2,col_pp1) = fun([211 211 211],tint_factor_2(2)); %light gray
% col_poi(3,col_pp1) = fun([105 105 105],tint_factor_2(3)); %dim gray
col_poi(1,col_pp1) = fun([211 211 211],tint_factor_2(1)); %light gray
col_poi(2,col_pp1) = fun([105 105 105],tint_factor_2(2)); %dim gray
col_poi(3,col_pp1) = fun([211 211 211],tint_factor_2(3)); %light gray
%col_poi(4,col_pp1) = fun([105 105 105],tint_factor_2(4)); %dim gray
%col_poi(5,col_pp1) = fun([105 105 105],tint_factor_2(5)); %dim gray
% col_poi(1,1:3) = fun([211 211 211]); %light gray
% col_poi(2,1:3) = fun([105 105 105]); %dim gray

%col_poi(1,col_pp2) = fun([0 0 0],tint_factor_2(1)); %black
col_poi(1,col_pp2) = fun([211 211 211],tint_factor_2(1)); %light gray
col_poi(2,col_pp2) = fun([105 105 105],tint_factor_2(2)); %dim gray
col_poi(3,col_pp2) = fun([211 211 211],tint_factor_2(3)); %light gray
%col_poi(4,col_pp2) = fun([105 105 105],tint_factor_2(4)); %dim gray
%col_poi(5,col_pp2) = fun([105 105 105],tint_factor_2(5)); %dim gray
% % col_poi(1,4:6) = fun([211 211 211]); %light gray
% col_poi(2,4:6) = fun([105 105 105]); %dim gray

%transparency of the face color
falpha = [1 1 1 1 1 ];
%transparency of the edge color
ealpha = [0.2 0.2 0.2 0.2 0.2];


%% prepare data for plotting
data = cell(1,numel(groups));
ave_data = NaN(1,numel(groups));
std_data = NaN(1,numel(groups));
sem_data = NaN(1,numel(groups));
er = cell(1,numel(groups));
X_rep = cell(1,numel(groups));
X_swarm = cell(1,numel(groups));
Y_swarm = cell(1,numel(groups));
max_val = NaN(1,numel(groups));
min_val = NaN(1,numel(groups));

% positions on the X axis
X = (1:numel(groups));

% data readout and stats calculation
for gi = 1:numel(groups)
    if plot_mode == 0
    
        data{1,gi} = eval(strcat(groups{gi},'.',parameter));
        data{1,gi} = data{1,gi}.*scale;
    elseif plot_mode == 1
       
        data_temp = eval(strcat(groups{gi},'.',parameter));
        
        if rheo_mode == 1
            rheo = eval(strcat(groups{gi},'.','Rheobase'));
        elseif rheo_mode == 2
            rheo = eval(strcat(groups{gi},'.','Rheobase_ind'));
        end
        
        data{1,gi} = get_step_values(data_temp, inc, rheo, rheo_mode, select_mode, stim).*scale;
        
    end
    
    
    ave_data(1,gi) = mean(data{1,gi},'omitnan');
    std_data(1,gi) = std(data{1,gi},'omitnan');
    sem_data(1,gi) = nansem(data{1,gi});
    
    max_val(1,gi) = max(data{1,gi});
    min_val(1,gi) = min(data{1,gi});
    
    X_rep{1,gi} = repmat(X(1,gi),1,numel(data{1,gi}));
    [X_swarm{1,gi}, Y_swarm{1,gi}] = swarmplot(X_rep{1,gi}+dp_offset,data{1,gi},0.05);
end

%% plotting
if numel(groups) <=3
    figure('Position',[500 100 350 500])
else
    figure('Position',[500 100 450 500])
    %figure('Position',[500 100 600 500]);
end

 
for pi = 1:numel(groups)
    
    if ~box_on

        if isnan(ave_data(1,pi))
            continue
        else
            bar(X(1,pi), ave_data(1,pi),... 
                'EdgeColor', col_bar{pi}, ...
                'EdgeAlpha', ealpha(pi),...
                'FaceColor', col_bar{pi}, ...
                'FaceAlpha', falpha(pi),...
                'LineWidth',2, 'BarWidth', 0.4,...
                'Clipping','off')
            hold on

            er{1,pi} = errorbar(X(1,pi), ave_data(1,pi), sem_data(1,pi));
            er{1,pi}.Color = 'k';
            er{1,pi}.LineWidth = 2;
            er{1,pi}.CapSize = 12;
        end
        
    else
            boxchart(X_rep{1,pi}(1,1:end)',data{1,pi}(1:end,1),...
                'BoxFaceColor',col_bar{pi},...
                'WhiskerLineColor',col_bar{pi},...
                'BoxFaceAlpha',ealpha(pi),...
                'LineWidth', 2,...
                'BoxWidth', 0.5,...
                'MarkerStyle','none',...
                'Notch',notch_on)
             
            hold on      
    end

    scatter(X_swarm{pi},Y_swarm{pi},...,
        'MarkerEdgeColor',col_poi(pi,4:6),...
        'MarkerFaceColor',col_poi(pi,1:3),...
        'LineWidth',1, ...
        'SizeData',30)
%             swarmchart(X_rep{1,pi}+dp_offset,data{1,pi},'XJitter','density',...
%             'XJitterWidth',0.2,'MarkerEdgeColor',col_poi(pi,4:6),...
%             'MarkerFaceColor',col_poi(pi,1:3),'LineWidth',1,'SizeData',20) 
end
      
%aesthetics of the plot
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;
ax.YLabel.String = strcat(parlabel,' (',unit,')');
ax.YAxis.Exponent = y_axis_expo;
ax.YDir = y_dir;

if strcmp(y_dir, y_dir_all{2})
    ax.YLim = [min(min_val)*1.15 y_lo_lim];
else
    ax.YLim = [-y_lo_lim max(max_val)*1.15];
end

%ax.YTick = [0 50 100];

ax.XTick = X;
%ax.XTick = [];
% ax.XTickLabel = [];
ax.XTickLabel = g_label;
ax.TickLabelInterpreter = 'none';
ax.XAxisLocation = 'origin';
% ax.YTick = [0 5 10 15];
% ax.YLim = [0 15];
%ax.XTickLabel.FontSize = 14;
%ax.XTickLabelRotation = 60;
%legend('hM4Di+CNO','Scaled','ctrl+CNO','FontSize',12,'Location','southeast')


hold off
box off

if plot_mode == 1 && select_mode == 1
     title(strcat('current injection: ',num2str(stim),'pA'))
end