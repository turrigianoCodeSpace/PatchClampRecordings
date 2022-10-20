function [vals] = get_step_values(data, inc, rheobase, rheo_mode, mode, stim)
% This function takes a processed dataset derived from f-I analysis and
% return the value(s) of the desired trace from all cells

% Inputs:
% data- the processed dataset (should not use raw data here)
% mode- 
    % mode=1: no normalization (starting from the 1st current step)
    % mode=2: normalized to the rheobase (starting from rheobase step)
% inc- current injection increment
% stim- stimulation amplitude (current injected, in pA)
% rheobase- rheobase values of the dataset (should also in the same mat
% file)
% rheo_mode -
     % 1: rheobase values in pA
     % 2: rheobase ind obtained from the fI data file

% Outputs:
% vals- values from the indicated traces (from all cells, stored in a
% column)

cell_num = size(data,2); % number of cells in this dataset

% output data
vals = NaN(cell_num,1);

rheo_ind = NaN(1,cell_num);
trace_ind = NaN(1,cell_num);

for ci = 1:cell_num
    
    
    if rheo_mode == 1
        rheo_ind(1,ci) = rheobase(ci,1)/inc;
    elseif rheo_mode == 2
        rheo_ind(1,ci) = rheobase(ci,1);
    end

    
    if mode == 1
        if stim == 0
            warning('Trace 0 does not exist!')
        else
            trace_ind(1,ci) = stim/inc;
        end
    elseif mode == 2
        trace_ind(1,ci) = rheo_ind(1,ci) + stim/inc;
    else
        warning('Unknown mode')
    end

    if isnan(trace_ind(1,ci)) %|| trace_ind(1,ci) > 20
        vals(ci,1) = NaN;
        warning(strcat('Cell', num2str(ci),' Trace', num2str(trace_ind(1,ci)),' does not exist!'))
    else

        vals(ci,1) = data(trace_ind(1,ci),ci);

    end
    
end
    
end


