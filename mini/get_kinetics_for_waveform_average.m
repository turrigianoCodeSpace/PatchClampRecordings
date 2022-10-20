%%%% This script calculates the kinetics and event charge of waveform averages 
%%%% of mini events obtained from a specific cell. 
%%%% WAVG and wavgs_raw yielded from MINIANALYSIS should be
%%%% loaded into the workspace before running.
%% change these for each run
% experiment type
% 1 for mEPSC, 2 for mIPSC
expt_type = 1;

%location of MINIANALYSIS mat files
fp_mini_analysis = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\culture_experiments\analyzed_mini_results\';

%rise time cutoff
if expt_type == 1
    rise = 'rise_1';
elseif expt_type == 2
    rise = 'mIPSC';
end
    
%Name of experiment(s) to be run
experiment = '220903';

%location where kinetic results are saved
fp_kinetics = ...,
    'C:\Users\schum\Google_Drive\Lab\Data_analysis\culture_experiments\waveform_average_kinetics\';

%%%%%plotting settings
figure_on = 0;

%%%%save results
save_results = 1;

%calc_mode = 1: calculate the average waveform across all traces of the
%same cell and then get kinetics of this wave. calc_mode = 2: calculate
%kinetics for each trace and average them across the same cell.
calc_mode = 1; 

%%

cd(strcat(fp_mini_analysis,rise))
load(strcat('MINIANALYSIS_',experiment,'.mat'))

%initiation
riseAve = cell(1,numel(WAVG));
decayAve = cell(1,numel(WAVG));
decayTau_1 = cell(1,numel(WAVG));
decayTot_1 = cell(1,numel(WAVG));
decayTau_2 = cell(1,numel(WAVG));
decayTot_2 = cell(1,numel(WAVG));
wavg_per_cell = cell(1,numel(WAVG));
wavg_charge = cell(1,numel(WAVG));
coefs_1 = cell(1,numel(WAVG));
gofs_1 = cell(1,numel(WAVG));
coefs_2 = cell(1,numel(WAVG));
gofs_2 = cell(1,numel(WAVG));
p1s = cell(1,numel(WAVG));
p2s = cell(1,numel(WAVG));
component_percentage = cell(1,numel(WAVG));

%Parameter settings
sample_rate = 5000;
rise_start_pc = 0.1;
rise_end_pc = 0.9;
decay_start_pc = 0.9;
num_pts_BL = 30;

if expt_type == 1
    num_pts_for_decayfit = 50;
elseif expt_type == 2
    num_pts_for_decayfit = 150;
end

%%%%Fitting parameters for decay phase
s1 = fitoptions('Method','NonlinearLeastSquares','Startpoint',[11 -0.4],...
    'Lower',[5 -3],'Upper',[120 -0.04]);
s2 = fitoptions('Method','NonlinearLeastSquares','Startpoint',[11 -0.4 5 -0.2],...
    'Lower',[5 -3 2 -3],'Upper',[120 -0.04 100 -0.01]);

f1 = fittype('a*exp(b*x)','options',s1);
f2 = fittype('a*exp(b*x)+c*exp(d*x)','options',s2);


for ei = 1:numel(wavgs_raw)
    for ci = 1:numel(wavgs_raw{1,ei})
        wavg_all_per_cell = [];
        
        if isempty(WAVG{1,ei}{1,ci})
            continue
        end
        
        if calc_mode == 1
            
            for ti = 1:numel(wavgs_raw{1,ei}{1,ci})
                if isempty(wavgs_raw{1,ei}{1,ci}{1,ti})
                    continue
                elseif isnan(wavgs_raw{1,ei}{1,ci}{1,ti})
                    continue
                else
                
                    for mi = 1:size(wavgs_raw{1,ei}{1,ci}{1,ti},1)
                        for ni = 1:size(wavgs_raw{1,ei}{1,ci}{1,ti},2)

                            if wavgs_raw{1,ei}{1,ci}{1,ti}(mi,ni) == 0
                                wavgs_raw{1,ei}{1,ci}{1,ti}(mi,ni) = NaN;
                            end
                        end
                    end

                    wavg_all_per_cell = cat(2,wavgs_raw{1,ei}{1,ci}{1,ti});
                end
            end

            for ii = 1:size(wavg_all_per_cell,1)
                wavg_per_cell{1,ei}(ii,ci) = mean(wavg_all_per_cell(ii,:),'omitnan');
            end
            
        elseif calc_mode ==2
            wavg_all_per_cell = WAVG{1,ei}{1,ci};
            
            for ii = 1:size(wavg_all_per_cell,1)
                wavg_per_cell{1,ei}(ii,ci) = mean(wavg_all_per_cell(ii,:),'omitnan');
            end
        end
        
        current_cell_average = wavg_per_cell{1,ei}(:,ci);
        if figure_on == 1
            figure('position',[56 200 1200 490])
            subplot(1,3,1:2)
            plot(0:0.2:0.2*(numel(current_cell_average)-1),current_cell_average)
            %plot(current_cell_average)
            title(strcat('Cell_',num2str(ci)),'Interpreter','none')
            ylabel('pA')
            xlabel('ms')
            hold on
        end


        %%%%calculate risetime: use interpolation to get a better resoultion
        
        [amplitude,pk_ind] = min(current_cell_average);
        
        if pk_ind > 100
            [amplitude,pk_ind] = min(current_cell_average(1:pk_ind-20));
        end
        
        rise_raw = current_cell_average(pk_ind-14:pk_ind); %use slope to find starting point
        slope = NaN(1,15-2);
        for si = 1:15-2 %calculate slope every 3 points
            slope(si) = (rise_raw(si+2)-rise_raw(si))/2;
        end
        
        
        
        interp_start = find(slope < -0.1,1,'first');
        
        interp_end = numel(rise_raw);
        
        interp_rise = interp1(interp_start:interp_end,rise_raw(interp_start:interp_end),interp_start:0.1:interp_end);
        
        rise_start_ind = find(interp_rise <= rise_start_pc*amplitude,1,'first');
        
        %recalibrate rise start point when the baseline is too tilted 
        if (numel(interp_rise) - rise_start_ind) > 60 
            [max_val,max_ind] = max(interp_rise(rise_start_ind : rise_start_ind+60));
            bl_med = median(interp_rise(rise_start_ind : rise_start_ind+60),'omitnan');
        
            if max_ind ~= 1
                f = fit((1:max_ind)',interp_rise(rise_start_ind:rise_start_ind+max_ind-1)','poly1');
                rise_slope = coeffvalues(f);

                if rise_slope(1) > 0.01
                   rise_start_ind = rise_start_ind + max_ind + ...,
                       find(interp_rise(rise_start_ind+max_ind: end)<=...,
                       bl_med-rise_start_pc*(bl_med-amplitude),1,'first');
                end
            end
        end
            
        rise_end_ind = find(interp_rise >= rise_end_pc*amplitude,1,'last');
        rise_time = (rise_end_ind-rise_start_ind)/(10*sample_rate);
        interp_time_interval = 1000/(sample_rate/0.1); %in ms
        
        rise_time_stamp = rise_start_ind*interp_time_interval : interp_time_interval : rise_end_ind*interp_time_interval;
        rise_time_stamp = rise_time_stamp + (pk_ind-14-1+interp_start-1)*0.2;
        if figure_on == 1
            plot(rise_time_stamp,interp_rise(rise_start_ind:rise_end_ind),'r','LineWidth',3)
        end
        
        riseAve{1,ei}(ci,1) = rise_time; %in s

        %%%%calculate decaytime
        decay_start_ind = find(current_cell_average(pk_ind:end) >= decay_start_pc*amplitude,1,'first')+pk_ind-1;
        decay_end_ind = decay_start_ind + num_pts_for_decayfit;

        decay_phase = current_cell_average(decay_start_ind:decay_end_ind);
        decayAve{1,ei}(ci,1) = numel(decay_phase)/sample_rate; %in s
        
        if figure_on == 1
            plot(decay_start_ind*0.2:0.2:decay_end_ind*0.2,decay_phase,'b','LineWidth',3)
        end

        %%%%Exponential fitting of decay phase
        %%%% the goodness of fit will be checked by a KS test, not R^2
        
        %first-order 

        [exp_fit_1,gof_1] = fit((0:0.2:0.2*num_pts_for_decayfit)',-decay_phase,f1);
        coefs_1{1,ei}(ci,1) = exp_fit_1.a;
        coefs_1{1,ei}(ci,2) = exp_fit_1.b;
        gofs_1{1,ei}(ci) = gof_1;
        preval_1 = feval(exp_fit_1,0:0.2:0.2*num_pts_for_decayfit);
        [h1,p1,ks2stat1] = kstest2(-decay_phase,preval_1);
        p1s{1,ei}(ci,1) = p1;
        
        %second-order 

        [exp_fit_2,gof_2] = fit((0:0.2:0.2*num_pts_for_decayfit)',-decay_phase,f2);
        coefs_2{1,ei}(ci,1) = exp_fit_2.a;
        coefs_2{1,ei}(ci,2) = exp_fit_2.b;
        coefs_2{1,ei}(ci,3) = exp_fit_2.c;
        coefs_2{1,ei}(ci,4) = exp_fit_2.d;
        gofs_2{1,ei}(ci) = gof_2;
        preval_2 = feval(exp_fit_2,0:0.2:0.2*num_pts_for_decayfit);
        [h2,p2,ks2stat2] = kstest2(-decay_phase,preval_2);
        p2s{1,ei}(ci,1) = p2;
        

        if figure_on == 1
            subplot(1,3,3)
            plot(0:0.2:0.2*num_pts_for_decayfit,-decay_phase)
            hold on
   
            plot(exp_fit_1,'r')
            plot(exp_fit_2,'c')
            title('decay phase fitting')
            legend('data','1st-order exp','2nd-order exp')
            ylabel('pA')
            xlabel('ms')
            hold off

        end

        decayTau_1{1,ei}(ci,1) = 1/(-exp_fit_1.b*1000); %in s
        decayTot_1{1,ei}(ci,1) = decayTau_1{1,ei}(ci,1)*log(exp_fit_1.a); %decay time calculated by exponential fit, in s
        
        decayTau_2{1,ei}(ci,1) = 1/(-exp_fit_2.b*1000); %fast 
        decayTau_2{1,ei}(ci,2) = 1/(-exp_fit_2.d*1000); %slow
        decayTot_2{1,ei}(ci,1) = decayTau_2{1,ei}(ci,1)*log(exp_fit_2.a); %fast
        decayTot_2{1,ei}(ci,2) = decayTau_2{1,ei}(ci,2)*log(exp_fit_2.c); %slow
        
        %for the 2nd exp fit, calculate the percentage of fast and slow
        %components at the peak (t=0)
        
        component_percentage{1,ei}(ci,1) = exp_fit_2.a/(exp_fit_2.a + exp_fit_2.c); %fast
        component_percentage{1,ei}(ci,2) = exp_fit_2.c/(exp_fit_2.a + exp_fit_2.c); %slow
        
        

        %%%%calculate charge
        wavg_charge{1,ei}(ci,1) = sum(current_cell_average((pk_ind-14-1+interp_start-1):decay_end_ind),'omitnan')...,
            /(-10^12*sample_rate); %in C
        if figure_on == 1
            subplot(1,3,1:2)
            for pi = (pk_ind-14-1+interp_start-1):decay_end_ind
                plot([pi*0.2 pi*0.2],[0 current_cell_average(pi)],'k-')
            end
        end

    end
    filename = strcat(experiment,'_calc_mode',num2str(calc_mode),'_',rise,'.mat');
end

%% save results
if save_results == 1
    cd(strcat(fp_kinetics,rise))
    save(filename,'wavg_per_cell','riseAve','decayTau_1','decayTot_1','decayTau_2','decayTot_2',...,
        'wavg_charge','coefs_1','gofs_1','coefs_2','gofs_2','p1s','p2s','component_percentage','cell_id')
end
