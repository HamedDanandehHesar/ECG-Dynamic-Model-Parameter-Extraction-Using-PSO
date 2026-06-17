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

% Nonlinear Phase using Dynamic Time warping
nonlinearphase = calculate_dtw_phase(ecg, qrs_positions);   % DTW-based nonlinear phase

%% -------- ECG mean and standard deviation in phase domain
[ECGsd,ECGmean,meanphase] = ecgsd_extractor_ver1(ecg,nonlinearphase,ecg_bins);

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
    dTheta = rem(nonlinearphase - Theta_i(j) + pi, 2*pi) - pi;

    % Add Gaussian contribution
    Z = Z + Alpha_i(j) .* exp(-dTheta.^2./(2*Beta_i(j).^2));
end

t = 1:length(ECG);

figure(3)
plot(t,ECG,'b', t,Z,'--r')
legend({'Original ECG', 'Synthetic ECG'})

OptimumParams = [Alpha_i Beta_i Theta_i];

savefile_params  = [file(1:end-4) '_params_nonlinear_Phase.mat']; 
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


function Phase = calculate_dtw_phase(ecg, r_peaks)

N = length(ecg);            % Total signal length
Phase = zeros(1,N);         % Phase array

nBeats = length(r_peaks);

%% -------- Extract beats centered on R-peaks
beats = {};                 % Cell array of beats
beat_index = {};            % Indices of each beat in original ECG
L = [];                     % Beat lengths

for i = 2:nBeats-1

    RR_prev = r_peaks(i) - r_peaks(i-1);   % Previous RR interval
    RR_next = r_peaks(i+1) - r_peaks(i);   % Next RR interval

    % Beat boundaries centered on R-peak
    start_idx = round(r_peaks(i) - RR_prev/2);
    end_idx   = round(r_peaks(i) + RR_next/2);

    start_idx = max(1,start_idx);
    end_idx   = min(N,end_idx);

    beat = ecg(start_idx:end_idx);

    beats{end+1} = beat;
    beat_index{end+1} = start_idx:end_idx;
    L(end+1) = length(beat);
end

%% -------- Build average template beat
Lref = round(mean(L));          % Reference length
template = zeros(1,Lref);

for i=1:length(beats)
    b = resample(beats{i},Lref,length(beats{i}));
    template = template + b;
end

template = template / length(beats);

%% -------- Define template phase (R at center)
phase_template = linspace(-pi,pi,Lref);

%% -------- DTW alignment and phase mapping
for i=1:length(beats)

    beat = beats{i};

    % DTW alignment indices
    [~,ix,iy] = dtw(beat,template);

    phi = zeros(1,length(beat));

    % Assign phase using warping path
    for k=1:length(ix)
        phi(ix(k)) = phase_template(iy(k));
    end

    % Interpolate missing phase samples
    phi = interp1(1:length(phi),phi,1:length(phi),'linear','extrap');

    Phase(beat_index{i}) = phi;
end

%% -------- Phase extrapolation before first beat
first_idx = beat_index{1}(1);
RR_mean = mean(diff(r_peaks));
stepTheta = 2*pi/RR_mean;

theta = Phase(first_idx);

for j = first_idx-1:-1:1
    theta = theta - stepTheta;
    theta = mod(theta+pi,2*pi)-pi;
    Phase(j) = theta;
end

%% -------- Phase extrapolation after last beat
last_idx = beat_index{end}(end);
theta = Phase(last_idx);

for j = last_idx+1:N
    theta = theta + stepTheta;
    theta = mod(theta+pi,2*pi)-pi;
    Phase(j) = theta;
end
end