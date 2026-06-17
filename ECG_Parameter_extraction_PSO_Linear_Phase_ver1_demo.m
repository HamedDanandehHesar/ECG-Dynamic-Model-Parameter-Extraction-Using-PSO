clear                     % Remove all variables from workspace
close all                 % Close all open figure windows

% Load ECG data from MAT file
[file, path] = uigetfile('*.mat','Select ECG mat file');
data = load(file);


fs = data.fs;                % Sampling frequency of ECG signal (Hz)

x = data.x;               % Extract stored signal matrix
ecg = x(1,:);             % Use first channel as ECG signal
length_sig = length(ecg); % Total number of ECG samples

ecg_bins = 250;           % Number of phase bins for ECG mean calculation

%% -------- R-peak detection using Pan–Tompkins algorithm
[qrs_positions] = pantompkins_qrs(ecg,fs);

figure(1),plot(ecg,'b'),hold on,plot(qrs_positions,ecg(qrs_positions),'*r'),hold off
legend({'ECG Signal','R Peaks'})
title(file)
%% -------- Phase calculation
% Linear phase based on RR intervals
[phase,~] = calculate_linear_phase_ver2(qrs_positions,length_sig,fs);

% Alternative options:
% phase = calculate_dtw_phase(ecg, qrs_positions);   % DTW-based nonlinear phase

%% -------- ECG mean and standard deviation in phase domain
[ECGsd,ECGmean,meanphase] = ecgsd_extractor_ver1(ecg,phase,ecg_bins);

%% -------- ECG parameter extraction using Gaussian mixture model

MaxNumGaussian = 25;   % Maximum number of Strongest Gaussian components

% ========================building new myfun based on L Gaussians
L_num_of_Gaussian_kernels = 50;
ecg_mean_temp = 0;
ai = [];
bi = [];
tetai  = [];
for i=1:L_num_of_Gaussian_kernels
% disp(num2str(i))
ecg_mean_temp1 = ECGmean - ecg_mean_temp;
lb = [-1.5*max(ecg_mean_temp1).*ones(1,1)   0.000001*ones(1,1)   (-pi+.014)*ones(1,1)  ];
ub = [(1.5*max(ecg_mean_temp1)).*ones(1,1)  5*ones(1,1)  (pi-.014)*ones(1,1)  ];  
myfun1 = @(params)  norm(ecg_mean_temp1'-sum((repmat(params(1:1),ecg_bins,1).*exp(-(rem(repmat(meanphase,1,1)'-repmat(params(3),ecg_bins,1)+pi,2*pi)-pi) .^2 ./ (2*(repmat(params(2),ecg_bins,1)) .^ 2))),2));


% options = optimoptions('particleswarm','SwarmSize',30,'HybridFcn',@fmincon,'MaxIter',1000);
options = optimoptions('particleswarm','SwarmSize',50,'MaxIter',100,'Display','off');

OptimumParams = particleswarm(myfun1,3*1,lb,ub,options);

% L = (length(OptimumParams)/3);

ai_1 = OptimumParams(1);
bi_1 = OptimumParams(2);
tetai_1 = OptimumParams(3);
ai = [ai ai_1];
bi = [bi bi_1];
tetai  = [tetai tetai_1];
dtetai_1 = rem(meanphase - tetai_1 + pi,2*pi)-pi;
ecg_mean_temp = ecg_mean_temp + ai_1 .* exp(-dtetai_1 .^2 ./ (2*bi_1 .^ 2));
figure(41),plot(ecg_mean_temp,'b'),hold on,plot(ECGmean,'r')
legend({'Synthetic ECG','ECG Mean'}),hold off
title([num2str(i) 'th' '  Gaussian'])
% pause(3)
end

[~,indx_strongest_peaks] = sort(abs(ai),'descend');

ai = ai(indx_strongest_peaks(1:MaxNumGaussian));
bi = bi(indx_strongest_peaks(1:MaxNumGaussian));
tetai = tetai(indx_strongest_peaks(1:MaxNumGaussian));
Synthetic_ECG_mean = 0;
for j=1:MaxNumGaussian
    dtetai = rem(meanphase - tetai(j) + pi,2*pi)-pi;

    Synthetic_ECG_mean = Synthetic_ECG_mean+ai(j) .* exp(-dtetai .^2 ./ (2*bi(j) .^ 2));
end

%% -------- Plot ECG mean vs synthetic ECG mean
figure(2);
plot(Synthetic_ECG_mean,'b'), hold on
plot(ECGmean,'--r')
iptsetpref('ImshowBorder','tight')
legend({'Synthetic ECG mean','Original ECG mean'})
title(['Sum of Top ' num2str(MaxNumGaussian) ' Gaussians'])
hold off

%% -------- Generate synthetic ECG signal using phase
ECG = x(1,:);

Alpha_i = ai;   % Gaussian amplitudes
Beta_i  = bi;   % Gaussian widths
Theta_i = tetai;   % Gaussian centers

Z = zeros(1,length(ECG));       % Synthetic ECG signal

for j = 1:length(Alpha_i)
    % Phase difference wrapped to [-pi, pi]
    dTheta = rem(phase - Theta_i(j) + pi, 2*pi) - pi;

    % Add Gaussian contribution
    Z = Z + Alpha_i(j) .* exp(-dTheta.^2./(2*Beta_i(j).^2));
end

t = 1:length(ECG);

figure(3)
plot(t,ECG,'b', t,Z,'--r')
legend({'Original ECG', 'Synthetic ECG'})

OptimumParams = [Alpha_i Beta_i Theta_i];

savefile_params  = [file(1:end-4) '_params_linear_Phase.mat']; 
save(savefile_params,'x',"OptimumParams","ECGmean","fs")



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions

function [ecgsd,ecg_mean,phase_mean] = ecgsd_extractor_ver1(ecg,phase,bins)

x1 = ecg;                        % ECG signal
meanPhase = zeros(1,bins);       % Mean phase per bin
ECGmean = zeros(1,bins);         % Mean ECG per bin
ECGsd = zeros(1,bins);           % ECG standard deviation per bin

% Handle wrap-around phase bin near -pi / +pi
I = find( phase >= (pi-pi/bins) | phase < (-pi+pi/bins) );

if(~isempty(I))
    meanPhase(1) = -pi;
    ECGmean(1) = mean(x1(I));
    ECGsd(1) = std(x1(I));
else
    ECGsd(1) = -1;               % Mark empty bins
end

% Loop over phase bins
for i = 1 : bins-1
    I = find( phase >= 2*pi*(i-0.5)/bins - pi & ...
              phase <  2*pi*(i+0.5)/bins - pi );

    if(~isempty(I))
        meanPhase(i+1) = mean(phase(I));
        ECGmean(i+1) = mean(x1(I));
        ECGsd(i+1) = std(x1(I));
    else
        ECGsd(i+1) = -1;
    end
end

% Interpolate missing bins
K = find(ECGsd==-1);

for i = 1:length(K)
    switch K(i)
        case 1
            meanPhase(1) = -pi;
            ECGmean(1) = ECGmean(2);
            ECGsd(1) = ECGsd(2);
        case bins
            meanPhase(bins) = pi;
            ECGmean(bins) = ECGmean(bins-1);
            ECGsd(bins) = ECGsd(bins-1);
        otherwise
            meanPhase(K(i)) = mean(meanPhase([K(i)-1 K(i)+1]));
            ECGmean(K(i))   = mean(ECGmean([K(i)-1 K(i)+1]));
            ECGsd(K(i))     = mean(ECGsd([K(i)-1 K(i)+1]));
    end
end

phase_mean = meanPhase;
ecg_mean   = ECGmean;
ecgsd      = ECGsd;

end



function [Phase,Omega] = calculate_linear_phase_ver2(locs,length_sig,fs)

% locs       : indices of detected R-peaks
% length_sig : total number of ECG samples
% fs         : sampling frequency

ind = locs(:)';                % Convert R‑peak indices to row vector

Phase = zeros(1,length_sig);   % Phase of each ECG sample
Omega = zeros(1,length_sig);   % Instantaneous angular frequency

RR = mean(diff(ind));          % Mean RR interval (samples)

%% -------- Phase before the first R‑peak

stepTheta = 2*pi/RR;           % Average phase increment per sample
omega_val = fs*stepTheta;      % Instantaneous angular frequency

theta = 0;                     % Initialize phase

for j = ind(1)-1:-1:1          % Move backward from first R‑peak
    theta = theta - stepTheta; % Decrease phase
    theta = mod(theta+pi,2*pi)-pi; % Wrap phase into [-pi , pi]

    Phase(j) = theta;          % Store phase
    Omega(j) = omega_val;      % Store frequency
end

%% -------- Phase between consecutive R‑peaks

for k = 1:length(ind)-1

    bins = ind(k+1)-ind(k);    % Number of samples between R-peaks

    stepTheta = 2*pi/bins;     % Phase increment so phase spans one cycle
    omega_val = fs*stepTheta;  % Corresponding angular frequency

    theta = 0;
    Phase(ind(k)) = 0;         % Define phase at R‑peak as zero

    for j = ind(k)+1 : ind(k+1)-1
        theta = theta + stepTheta; % Linear phase progression
        if theta>pi
            theta = -pi;
        end
        Phase(j) = theta;
        Omega(j) = omega_val;
    end

    Phase(ind(k+1)) = 0;       % Next R‑peak also set to zero phase
end

%% -------- Phase after the last R‑peak

stepTheta = 2*pi/RR;           % Use mean RR again
omega_val = fs*stepTheta;

theta = 0;

for j = ind(end)+1:length_sig
    theta = theta + stepTheta; % Continue phase linearly
    theta = mod(theta+pi,2*pi)-pi;

    Phase(j) = theta;
    Omega(j) = omega_val;
end

end