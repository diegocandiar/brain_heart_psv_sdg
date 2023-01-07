function struct_sdg = compute_psv_sdg(data_eeg, pks_indx)


%% get heartbeats timings
%data
time = data_eeg.time{1};
t_heartbeats = time(pks_indx);
IBI = diff(t_heartbeats);
t_IBI = time(pks_indx(1:length(IBI)));

%% Poincare
wind = 15;
[CSI_out, CVI_out, t_out] = compute_CSI_CVI(IBI, t_IBI, wind);

%% frequency analysis
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.pad          = 'nextpow2';
cfg.foi          = 0:0.5:45;                         
cfg.t_ftimwin    = ones(length(cfg.foi),1).*2;  
cfg.toi          = '50%';  % 2s 50p,                  
freq = ft_freqanalysis(cfg, data_eeg);

%% remove nans
time2 = freq.time;
time3 = time2(2:length(time2)-2);
t1 = max([time3(1) t_out(1)]);
t2 = min([time3(end) t_out(end)]);
time3 = time3(time3 >=t1 & time3 <=t2);

cfg = [];
%cfg.frequency = [1 45];
cfg.latency = [time3(1) time3(end)];
freq = ft_selectdata(cfg, freq);

CSI = interp1(t_out,CSI_out,time3,'spline');
CVI = interp1(t_out,CVI_out,time3,'spline');

%% freq bands
freq_epoch = freq; freq_epoch = rmfield(freq_epoch,{'powspctrm','freq'});
freq_delta = freq_epoch;
freq_theta = freq_epoch;
freq_alpha = freq_epoch;
freq_beta = freq_epoch;
freq_gamma = freq_epoch;

f1a = 1; f1b = 4; f1 = freq.freq >= f1a & freq.freq <= f1b;
f2a = 4; f2b = 8; f2 = freq.freq >= f2a & freq.freq <= f2b;
f3a = 8; f3b = 12; f3 = freq.freq >= f3a & freq.freq <= f3b;
f4a = 12; f4b = 30; f4 = freq.freq >= f4a & freq.freq <= f4b;
f5a = 30; f5b = 45; f5 = freq.freq >= f5a & freq.freq <= f5b;

freq_delta.trial{1} = squeeze(trapz(freq.powspctrm(:,f1,:),2));  
freq_theta.trial{1} = squeeze(trapz(freq.powspctrm(:,f2,:),2));    
freq_alpha.trial{1} = squeeze(trapz(freq.powspctrm(:,f3,:),2));  
freq_beta.trial{1} = squeeze(trapz(freq.powspctrm(:,f4,:),2));   
freq_gamma.trial{1} = squeeze(trapz(freq.powspctrm(:,f5,:),2));
freq_delta.dimord = 'chan_time';
freq_theta.dimord = 'chan_time';
freq_alpha.dimord = 'chan_time';
freq_beta.dimord = 'chan_time';
freq_gamma.dimord = 'chan_time';


%% BHI
win_RR = 15; % seconds
time = time3;
timea = time(floor(win_RR/2) +1 : end-ceil(win_RR/2));
timed = time(win_RR+1 : end-win_RR) + floor(win_RR/2);


%% SDG model CSI CVI 
Fs = 1;

% delta
EEG_comp = freq_delta.trial{1}; 
[bhi_CSI_delta, bhi_CVI_delta, bhi_delta_CSI, bhi_delta_CVI, ~, ~] = model_psv_sdg(EEG_comp, IBI, t_IBI, CSI, CVI, Fs, time, wind); 

%theta
EEG_comp = freq_theta.trial{1}; 
[bhi_CSI_theta, bhi_CVI_theta, bhi_theta_CSI, bhi_theta_CVI, ~, ~] = model_psv_sdg(EEG_comp, IBI, t_IBI, CSI, CVI, Fs, time, wind); 

%alpha
EEG_comp = freq_alpha.trial{1}; 
[bhi_CSI_alpha, bhi_CVI_alpha, bhi_alpha_CSI, bhi_alpha_CVI, ~, ~] = model_psv_sdg(EEG_comp, IBI, t_IBI, CSI, CVI, Fs, time, wind); 

%beta
EEG_comp = freq_beta.trial{1}; 
[bhi_CSI_beta, bhi_CVI_beta, bhi_beta_CSI, bhi_beta_CVI, ~, ~] = model_psv_sdg(EEG_comp, IBI, t_IBI, CSI, CVI, Fs, time, wind); 

%gamma
EEG_comp = freq_gamma.trial{1}; 
[bhi_CSI_gamma, bhi_CVI_gamma, bhi_gamma_CSI, bhi_gamma_CVI, ~, ~] = model_psv_sdg(EEG_comp, IBI, t_IBI, CSI, CVI, Fs, time, wind); 

struct_sdg = struct;
struct_sdg.time = time;
struct_sdg.timea = timea;
struct_sdg.timed = timed;
struct_sdg.CSI = CSI;
struct_sdg.CVI = CVI;
struct_sdg.freq_delta = freq_delta;
struct_sdg.freq_theta = freq_theta;
struct_sdg.freq_alpha = freq_alpha;
struct_sdg.freq_beta = freq_beta;
struct_sdg.freq_gamma = freq_gamma;
struct_sdg.bhi_delta_CSI = bhi_delta_CSI ;
struct_sdg.bhi_delta_CVI = bhi_delta_CVI;
struct_sdg.bhi_theta_CSI = bhi_theta_CSI;
struct_sdg.bhi_theta_CVI = bhi_theta_CVI;
struct_sdg.bhi_alpha_CSI = bhi_alpha_CSI;
struct_sdg.bhi_alpha_CVI = bhi_alpha_CVI;
struct_sdg.bhi_beta_CSI = bhi_beta_CSI;
struct_sdg.bhi_beta_CVI = bhi_beta_CVI;
struct_sdg.bhi_gamma_CSI = bhi_gamma_CSI ;
struct_sdg.bhi_gamma_CVI = bhi_gamma_CVI;
struct_sdg.bhi_CSI_delta = bhi_CSI_delta;
struct_sdg.bhi_CVI_delta = bhi_CVI_delta;
struct_sdg.bhi_CSI_theta = bhi_CSI_theta;
struct_sdg.bhi_CVI_theta = bhi_CVI_theta;
struct_sdg.bhi_CSI_alpha = bhi_CSI_alpha;
struct_sdg.bhi_CVI_alpha = bhi_CVI_alpha;
struct_sdg.bhi_CSI_beta = bhi_CSI_beta ;
struct_sdg.bhi_CVI_beta = bhi_CVI_beta;
struct_sdg.bhi_CSI_gamma = bhi_CSI_gamma;
struct_sdg.bhi_CVI_gamma = bhi_CVI_gamma;

end

