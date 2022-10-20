function [sem] = nansem(data)
% This function calculate standard error of a data set after removing NaN
% data- input data set
% sem- standard error from the mean

sem = nanstd(data)/sqrt(count_non_nan(data));
end


