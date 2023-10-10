function [swarm_X,swarm_Y] = swarmplot(X, Y, spread_width)
% This function spread indivual data points using the beeswarm alogrithm.
% Each point will have an offset scaled to their magnitudes and spread as
% in gaussian distribution.
% Input: X- the x-values of the data (since these data are usually
%           categorical, X for all data points should be the same.)
%        Y- the y-values of the data
%        spread_width- the horizontal distance onto which the data points
%                      are spread.
%        bin_width- the width of each bin that the data are divided into.
% Output: swarm_X- the X values for swarmplot
%         swarm_Y- the Y values for swarmplot

if ~isvector(X) || ~isvector(Y)
    disp('X and Y must be vectors!')
else
    val_min = min(Y);
    val_max = max(Y);
    
    Y = sort(Y);
    bin_width_min = Y(find(Y>prctile(Y,10),1,'first'));
    bin_width_max = Y(find(Y<prctile(Y,90),1,'last'));
    
    % for data that cluster within a narrow interval, take percentile
    % rather than the acutal value
    if bin_width_max <= bin_width_min
        bin_width_max = prctile(Y,90);
        bin_width_min = prctile(Y,10);
    end
    
    bin_width = (bin_width_max-bin_width_min)/8;
    spread_rows = round((val_max - val_min)/bin_width);
    bin_idx = 0;
    swarm_x_temp = NaN(numel(Y),spread_rows+1);
    swarm_y_temp = NaN(numel(Y),spread_rows+1);
    
    for bi = val_min : bin_width : val_max
        bin_idx = bin_idx + 1;
        num_bin = 0;
        current_y = NaN(1,numel(Y));
        
        for yi = 1:numel(Y)
            if Y(yi) >= bi && Y(yi) < bi+bin_width
                num_bin = num_bin + 1;
                current_y(num_bin) = Y(yi);
            end

        end
        
        %if more than 1 values exist in this bin, spread.
        %current_y = nonzeros(current_y);
        current_y_temp = NaN(1,num_bin);
        current_x = NaN(1,num_bin);    
        if num_bin > 1
            center_x = (num_bin+1)/2;
                   
            center_y = max(current_y);
            
            current_y_temp(round(center_x)) = center_y;


            for xi = 1:num_bin
                if xi <= center_x
                    current_x(xi) = X(1)-(center_x - xi)*spread_width;
                else
                    current_x(xi) = X(1)+(xi - center_x)*spread_width;
                end
                
                swarm_x_temp(1:numel(current_x),bin_idx) = current_x;
            end
            
            %sort y in gaussian distribution
            for yii = 1:2:(num_bin-1)
                 current_y_temp((yii+1)/2) = current_y(yii);
                 if num_bin > 2
                     current_y_temp(num_bin-(yii+1)/2+1) = current_y(yii+1);
                 end
            end
            
            swarm_y_temp(1:numel(current_y_temp),bin_idx) = current_y_temp;    
            
        elseif num_bin == 1
            swarm_x_temp(1,bin_idx) = X(1);
            swarm_y_temp(1,bin_idx) = current_y(1);
        end
    end
   
    swarm_X = swarm_x_temp(~isnan(swarm_x_temp));
    swarm_Y = swarm_y_temp(~isnan(swarm_y_temp));
end
