function [CSI_out, CVI_out, t_out] = compute_CSI_CVI(RR, t_RR, wind)

Fs = 4;
time = t_RR(1) : 1 / Fs : t_RR(end);

%% first poincare plot

sd=diff(RR); 
SD01 = sqrt(0.5*std(sd)^2);
SD02 = sqrt(2*(std(RR)^2)-(0.5*std(sd)^2));

%% time varying SD
t1 = time(1);
t2 = t1 +  wind;
ixs = find(t_RR > t2);
nt = length(ixs)-1;

SD1 = zeros(1,nt);
SD2 = zeros(1,nt);
t_C = zeros(1,nt);

for k = 1 : nt
    i = ixs(k); 

    t2 = t_RR(i);
    t1 = t_RR(i)-wind;
    ix = find(t_RR >= t1 & t_RR<= t2);
    
    sd=diff(RR(ix)); 
    SD1(k) = sqrt(0.5*std(sd)^2);
    SD2(k) = sqrt(2*(std(RR(ix))^2)-(0.5*std(sd)^2));

    t_C(k) = t2;

end

SD1 = SD1 - mean(SD1) + SD01;
SD2 = SD2 - mean(SD2) + SD02;

% CVI = SD1.*SD2 * 100;
% CSI = SD2./SD1;


CVI = SD1 * 10;
CSI = SD2 * 10;

%%
t_out = t_C(1) : 1 / Fs : t_C(end);
CVI_out = interp1(t_C, CVI, t_out, 'Spline');
CSI_out = interp1(t_C, CSI, t_out, 'Spline');

end