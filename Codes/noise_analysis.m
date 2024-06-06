clear all
close all
clc

load("PID.mat")

time = MLS2EMExpData.time;
position = MLS2EMExpData.signals(1).values;
velocity = MLS2EMExpData.signals(2).values;
current = MLS2EMExpData.signals(3).values(:,1);

figure(1)
subplot(3,1,1); 
hold on
grid on
plot(time,position)
title('Position');
xlabel('Time');
ylabel('Amplitude');
hold off

subplot(3,1,2); 
hold on
grid on
plot(time,velocity)
axis([0 10 -1 1])
hold off

subplot(3,1,3); 
hold on
grid on
plot(time,current)
hold off
% I remove the first X s because of the raise time
% X = 1750;
% X=2000;
X = 2900;
% X = 2500;
pos = position(X:end);
vel = velocity(X:end);
curr = current(X:end);
t = time(X:end);

%% Calculate Statistical Measures:
%mean, variance, and standard deviation
mean_pos = mean(pos);
var_pos = var(pos);
std_dev_pos = std(pos);

mean_vel = mean(vel);
var_vel = var(vel);
std_dev_vel = std(vel);

mean_curr = mean(curr);
var_curr = var(curr);
std_dev_curr = std(curr);

%% compute the noise
noise_pos = pos - mean_pos;
noise_vel = vel - mean_vel;
noise_curr = curr - mean_curr;

% % plot of the noise
% figure(3)
% subplot(3,1,1); 
% hold on
% grid on
% plot(t,noise_pos)
% title('Position Noise');
% xlabel('Time');
% ylabel('Amplitude');
% hold off
% 
% subplot(3,1,2); 
% hold on
% grid on
% plot(t,noise_vel)
% title('Velocity Noise');
% xlabel('Time');
% ylabel('Amplitude');
% hold off
% 
% subplot(3,1,3); 
% hold on
% grid on
% plot(t,noise_curr)
% title('Current Noise');
% xlabel('Time');
% ylabel('Amplitude');
% hold off


% Standard deviations of each process, specified as a vector of length n with the standard deviations of each process.
Exp_sigma = [std_dev_pos std_dev_vel std_dev_curr];


corr_pv = corrcoef(noise_pos,noise_vel);
corr_pc = corrcoef(noise_pos,noise_curr);
corr_vc = corrcoef(noise_vel,noise_curr);

%Correlation matrix, specified as an n-by-n correlation coefficient matrix. 
ExpCorr = [1 corr_pv(1,2) corr_pc(1,2);
           corr_pv(1,2)  1   corr_vc(1,2);
           corr_pc(1,2)  corr_vc(1,2) 1];
ExpCovariance = corr2cov(Exp_sigma,ExpCorr);

%% outliners

% Calculate z-scores
z_scores = (noise_pos - mean(noise_pos)) / std_dev_pos;

% Define threshold
threshold = 3;

% Identify and remove outliers
outliers = abs(z_scores) > threshold;
data_cleaned = noise_pos(~outliers);

% Calculate quartiles and IQR
Q1 = quantile(noise_pos, 0.25);
Q3 = quantile(noise_pos, 0.75);
IQR = Q3 - Q1;

% Define bounds for outliers
lower_bound = Q1 - 1.5 * IQR;
upper_bound = Q3 + 1.5 * IQR;

% Identify outliers
outliers = (noise_pos < lower_bound) | (noise_pos > upper_bound);

%% Calculate NEW Statistical Measures:

% Normalize statistics
normalized_mean = mean_pos; % Mean doesn't need normalization
normalized_std = std_dev_pos / mean_pos; % Normalize standard deviation by dividing by mean
normalized_var = var_pos / (mean_pos^2); % Normalize variance by dividing by square of mean

%tutte le misure vengono di  nuovo uguali
mean_noise_pos = mean(noise_pos);
mean_noise_vel = mean(noise_vel);
mean_noise_curr = mean(noise_curr);


%% Plot histogram
figure(21)
histogram(position, 'Normalization', 'pdf',NumBins=150)
hold on
y = 0.009:0.0000001:0.011;
mu = mean_pos;
sigma = std_dev_pos;
f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(y,f,'LineWidth',1.5)

%% Frequency Domain Analysis: 
Fs = 1000; % Sampling frequency (Hz), adjust according to your simulation
N = length(noise_pos); % Length of the signal
f = Fs*(0:(N/2))/N; % Frequency vector

Y = fft(noise_pos);
P = abs(Y/N);
P = P(1:N/2+1);

figure(31)
plot(f, P);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Single-Sided Magnitude Spectrum of Noise Signal');


%% Filtering: 

freq_taglio = 60; % [Hz]
Ts = 0.001;
n= [1-2.72^(-Ts*2*pi*freq_taglio) 0];
d = [1 -2.72^(-Ts*2*pi*freq_taglio)];
H = tf(n,d,Ts);

a = freq_taglio*2*3.14*Ts;
na = [a 0];
da = [1 a-1];
Ha = tf(na,da,Ts);

nc = [1];
dc = [1/(2*pi*freq_taglio) 1];
Hc = tf(nc,dc);
Hc_d = c2d(Hc,Ts);


figure(35)
bode(H)
hold on
bode(Ha)
bode(Hc_d)
legend('IIR', 'EWMA','Discretized Low Pass')

