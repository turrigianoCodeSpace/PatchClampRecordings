function [Rin, Cm, mTau, fit_vals, trace_time, exp_fit] = get_PP_I_clamp(data, step_start, pulse, I_step, samprate, figure_on)

%Calculates input resistance (Rin) and membrane capacitance (Cm) for traces
%obtained under current clamp

%Inputs:
%   data- the region from the raw trace selected for analysis
%   step start- this is where the sealtest starts, should be in seconds
%   pulse- duration of test pulse (in seconds)
%   I_step- test pulse amplitude (in A) 
%   samprate- sampling rate, in Hz
%   figures_on = 1 for plotting

%Outputs:
%   Rin- Input resistance (in MOhm)
%   Cm- Membrane capacitance (in pF)
%   mTau- Membrane time constant (in ms)
%   fit_vals- regions for 1st-exp fitting
%   trace_time- time points for 1st-exp fitting
%   exp_fit- fitting results


data_V = data * 10^-3; % mV to V
dV_sec = diff(data_V).*(samprate/1000); %mv/sec


step_start_ind = step_start*samprate;

bl_vals = data_V(step_start_ind-2000 : step_start_ind); %baseline region
bl_vals = filloutliers(bl_vals, 'spline');%, 'movmedian', 20);

%bl_std = std(bl_vals, 'omitnan');
bl_dv = dV_sec(step_start_ind-2000 : step_start_ind);
%bl_dv =  filloutliers(bl_dv, 'spline');
if ~isempty(find(bl_dv > 0.01, 1, 'first'))
    noi_ind = find(bl_dv > 0.01, 1, 'first');
    bl_vals(noi_ind:noi_ind+50) = NaN;
end
% for bi = 1:numel(bl_dv)
%     if bl_dv(bi) 
%         bl_vals(bi) = NaN;
%     end
% end

trace_end_ind = step_start_ind + pulse*samprate;
trace_vals = data_V(step_start_ind : trace_end_ind); %regions selected for fitting

% figure();
% plot(trace_vals);

ss_vals = trace_vals(end-2000 : end);
ss_vals = filloutliers(ss_vals, 'spline', 'movmedian', 20);
% ss_std = std(ss_vals, 'omitnan');
% for si = 1:numel(ss_vals)
%     if ss_vals(si) < mean(ss_vals,'omitnan')-3*ss_std || ss_vals(si) > mean(ss_vals,'omitnan')+3*ss_std
%         ss_vals(si) = NaN;
%     end
% end


% Find Rin
V_ss = median(ss_vals, 'omitnan');
V_bl = median(bl_vals, 'omitnan');

delta_V = V_ss - V_bl;
Rin = delta_V/I_step*10^-6; %in MOhm

% Find membrane constant tau (in seconds)
if isempty(find(trace_vals(1 : 800) >= median(bl_vals, 'omitnan')+0.1*delta_V, 1, 'last'))
    fit_st = 1;
else
    fit_st = find(trace_vals(1 : 800) >= median(bl_vals, 'omitnan')+0.1*delta_V, 1, 'last'); % 10% drop as the stating point
end

[~,fit_end] = min(trace_vals);
trace_vals_adj = trace_vals(fit_st:fit_end);

% trace_vals_abs = abs(trace_vals);

%%%single exponential fitting
%Normalization of the fitting region
trace_vals_norm = abs(trace_vals_adj - V_bl);
fit_st_norm = find(trace_vals_norm > trace_vals_norm(end)*0.05, 1, 'first');

% if find(trace_vals_norm > trace_vals_norm(end)*0.95, 1, 'first') < 150
%     fit_end_norm = find(trace_vals_norm > trace_vals_norm(end)*0.95, 1, 'last')+100;
% else
    fit_end_norm = find(trace_vals_norm > trace_vals_norm(end)*0.95, 1, 'last');
%end

fit_vals = trace_vals_norm(fit_st_norm:fit_end_norm);
%fit_vals_test = fit_vals - 10*10^-4;

trace_time = 0:1/samprate:(numel(fit_vals)-1)/samprate;

tau_est = (find(fit_vals>fit_vals(end)*0.63, 1,'first'))/samprate; 

s = fitoptions('Method', 'NonlinearLeastSquares', ...
    'StartPoint', [fit_vals(end), tau_est*0.5]);%, ...
%      'Lower', [mean(fit_vals(end-1:end)*.6), tau_est*0.001],...
%      'Upper', [fit_vals(end)*1.2, tau_est*5]);
f = fittype('a*(1-exp(-x/b))','options',s);

% s2 = fitoptions('Method', 'NonlinearLeastSquares', ...
%     'StartPoint', [fit_vals(end), tau_est*0.5, fit_vals(end-7), tau_est*1.5],...
%      'Lower', [mean(fit_vals(end-1:end)*.6), tau_est*0.001, mean(fit_vals(end-5:end))*0.1, tau_est*0.01],...
%      'Upper', [fit_vals(end)*1.2, tau_est*5, fit_vals(end)*1.25, tau_est*15]);
% f2 = fittype('a*(1-exp(-x/b))+c*(1-exp(-x/d))','options',s2);

[exp_fit,gof] = fit(trace_time',fit_vals,f);
%[exp_fit_2,gof_2] = fit(trace_time',fit_vals,f2);

mTau = exp_fit.b; %in s
Cm = exp_fit.b/Rin*10^6; %in pF
% Cm_2 = exp_fit_2.d/Rin*10^6; %in pF
% Cp = exp_fit_2.b/500*10^6; %in pF (hypothetically pipette capacitance)

%single-exponential fit estimation
fit_est_vals = exp_fit.a*(1-exp(-trace_time/exp_fit.b));

% %second-exponential fir estimation
% fit_est_vals_exp2 = exp_fit_2.a*(1-exp(-trace_time/exp_fit_2.b))+...,
%                     exp_fit_2.c*(1-exp(-trace_time/exp_fit_2.d));
% fit_est_vals_exp2_fast  = exp_fit_2.a*(1-exp(-trace_time/exp_fit_2.b));
% fit_est_vals_exp2_slow = exp_fit_2.c*(1-exp(-trace_time/exp_fit_2.d));

% plot fitted curves
if figure_on
    figure()
    hold on
    plot(trace_time', fit_vals, 'Color','#377899', 'LineWidth',2)
    plot(trace_time', fit_est_vals', 'Color','#E6818B', 'LineWidth',2)
    %plot(trace_time', fit_est_vals_exp2', 'Color','#E6892E', 'LineWidth',2)
    
    legend('data','single-exp','Location','southeast')
    xlabel('time (s)')
    ylabel('membrane potential (V)')
    hold off
    
    
%     figure()
%     hold on
%     plot(trace_time', fit_vals, 'Color','#377899', 'LineWidth',2)
%     plot(trace_time', fit_est_vals_exp2_fast', 'Color','#7EE66A', 'LineWidth',2)
%     plot(trace_time', fit_est_vals_exp2_slow', 'Color','#E675C8', 'LineWidth',2)
%     
%     xlabel('time (s)')
%     ylabel('membrane potential (V)')    
%     legend('data','fast','slow')
%     hold off
%     
end
       