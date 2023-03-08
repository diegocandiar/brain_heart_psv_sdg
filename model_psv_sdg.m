function [CSI2B, CVI2B, B2CSI, B2CVI, tH2B, tB2H] = model_psv_sdg(EEG_comp, IBI, t_IBI, CSI, CVI, Fs, time, wind)
% The function model_psv_sdg computes the coupling coefficients between EEG
% in a defined frequency band with the spectral estimators of symathovagal activity
% in both directions.
% The inputs are:
% -EEG_comp: Matrix ChannelxTime containing the time-varying power
% -IBI: NON interpolated interbeat intervals in seconds
% -t_IBI: time of each IBI in seconds
% -Fs: EEG sample frequency
% -time: the time array for EEG data
% -wind: time window in seconds in which the model is estimated (wind = 15)
% EEG=sqrt(EEG) to ease model convergence, else the model use of real values
%
% The outputs are:
% -CS2B: time-varying coupling coefficients from CSI to EEG band
% -CV2B: time-varying coupling coefficients from CVI to EEG band
% -B2CS: time-varying coupling coefficients from EEG band to CSI
% -B2CV: time-varying coupling coefficients from EEG band to CVI
% -tH2B: time array of the coefficients from heart to brain
% -tB2H: time array of the coefficients from brain to heart

% Author: Diego Candia-Rivera (diego.candia.r@ug.uchile.cl)

%% Check if one-channel EEG data has inverted dimension

[Nch,Nt] = size(EEG_comp);
if Nt==1
    EEG_comp = EEG_comp';
    [Nch,Nt] = size(EEG_comp);
end

%% HRV model based on Poincare plot
f_lf = 0.1;
f_hf = 0.25;

w_lf = 2*pi*f_lf;
w_hf = 2*pi*f_hf;

ss = wind*Fs;
sc = 1;
nt = ceil((length(time)-ss)/sc);

Cs = zeros(1,nt);
Cp = zeros(1,nt);
TM = zeros(1,nt);


for i = 1:nt
    ix1 = (i-1)*sc + 1;
    ix2 = ix1 + ss - 1;
    ixm = floor(mean(ix1:ix2));   
%     ixm = ix2;
    t1 = time(ix1);
    t2 = time(ix2);
    ix = find(t_IBI >= t1 & t_IBI<= t2);

    
    mu_ibi = mean(IBI(ix));
    mu_hr = 1/mu_ibi;  
    
    G = sin(w_hf/(2*mu_hr))-sin(w_lf/(2*mu_hr)); 
    
    M_11 = sin(w_hf/(2*mu_hr))*w_lf*mu_hr/(sin(w_lf/(2*mu_hr))*4);
    M_12 = -sqrt(2)*w_lf*mu_hr/(8*sin(w_lf/(2*mu_hr)));
    M_21 = -sin(w_lf/(2*mu_hr))*w_hf*mu_hr/(sin(w_hf/(2*mu_hr))*4);
    M_22 = sqrt(2)*w_hf*mu_hr/(8*sin(w_hf/(2*mu_hr)));
    M = [M_11, M_12; M_21, M_22];
    L = max(IBI(ix))-min(IBI(ix));         
    W = sqrt(2)*max(abs(IBI(ix(2:end))-IBI(ix(1:end-1))));
    C = 1/G*M*[L; W];
    Cs(i) = C(1);   Cp(i) = C(2);
    TM(i) = time(ixm);
end


%% interpolation (edges are extended to avoid extrapolation)
t1 = max([time(1) TM(1)]);
t2 = min([time(end) TM(end)]);
time2 = time(time >=t1 & time <=t2);

Cs = interp1(TM, Cs , time2, 'spline');
Cp = interp1(TM, Cp , time2, 'spline');

CSI = interp1(time, CSI , time2, 'spline');
CVI = interp1(time, CVI , time2, 'spline');

EEG_old = EEG_comp;
[Nch, Nt] = size(EEG_old);
clear EEG_comp

for i = 1 : Nch
    EEG_comp(i,:) = interp1(time, EEG_old(i,:), time2, 'spline');
end

time = time2;

% Optional
EEG_comp = clean_artif(EEG_comp);
% EEG_comp = sqrt(EEG_comp);
% Cs = Cs / std(Cs);
% Cp = Cp / std(Cp);

%% RUN MODEL

parfor ch = 1:Nch % Consider use parallel computing here
    [CSI2B(ch,:), CVI2B(ch,:), B2CSI(ch,:), B2CVI(ch,:)] = SDG(EEG_comp(ch,:), CSI, CVI, Cs, Cp, wind*Fs);
end

%% TIME VECTOR DEFINITION
% at the beginning of the window
tH2B = time(1 : end-wind*Fs);
tB2H = time(wind*Fs+1 : end-wind*Fs);

% in the middle of the window
% tH2B = time(floor(wind/2) +1 : end-ceil(wind/2));
% tB2H = time(wind+1 : end-wind) + wind/2;

end

function [CSI_to_EEG, CVI_to_EEG, EEG_to_CSI, EEG_to_CVI] = SDG(EEG_ch, HRV_CSI, HRV_CVI, Cs_i, Cp_i, window)

    Nt = length(EEG_ch);
    %% First time window is calculated separately
    for i = 1 : window
        arx_data = iddata(EEG_ch(i:i+window)', HRV_CSI(i:i+window)',1); 
        model_eegP = arx(arx_data,[1 1 1]);                                                 
        CSI_to_EEG(i) = model_eegP.B(2);

        arx_data = iddata(EEG_ch(i:i+window)', HRV_CVI(i:i+window)',1); 
        model_eegP = arx(arx_data,[1 1 1]);                                                 
        CVI_to_EEG(i) = model_eegP.B(2);
        
        pow_eeg(1,i) = mean(EEG_ch(i:i+window));                                 
    end
    
    
    for i = window+1:min([length(Cp_i),Nt-window, length(HRV_CSI)-window])

        %% Heart to brain estimation
        arx_data = iddata(EEG_ch(i:i+window)', HRV_CSI(i:i+window)',1); 
        model_eegP = arx(arx_data,[1 1 1]); 
        CSI_to_EEG(i) = model_eegP.B(2); 

        arx_data = iddata(EEG_ch(i:i+window)', HRV_CVI(i:i+window)',1); 
        model_eegP = arx(arx_data,[1 1 1]); 
        CVI_to_EEG(i) = model_eegP.B(2); 
        
        pow_eeg(1,i) = mean(EEG_ch(i:i+window));

        %% Brain to heart estimation
        if i-window <= length(Cp_i)-window-1
            EEG_to_CVI(i-window) = mean((Cp_i(i-window:i))./pow_eeg(i-window:i));
            EEG_to_CSI(i-window) = mean((Cs_i(i-window:i))./pow_eeg(i-window:i));
        else
            EEG_to_CVI(i-window) = EEG_to_CVI(i-window-1);
            EEG_to_CSI(i-window) = EEG_to_CSI(i-window-1);
        end
    end

end