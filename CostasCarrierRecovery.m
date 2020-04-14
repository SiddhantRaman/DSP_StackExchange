%% Costas Loop Carrier Phase Recovery : Named after J.P. Costas
clc
close all
clear all

%% Generate received signal $s[n] = cos(2\pi f_c k T_s + \phi)$
% $\phi$ is unknown phase of the carrier, which needs to be estimated

% This simulation assumes $\phi = 0.2$, $f_c=1000Hz$ and $f_s = 5000Hz$
fc = 50;                      % Carrier Frequency
Ts = 1/5000;                    % Sampling Period
n = 0:2000000-1;                % n
phi = 0.2;
s_n = 10 + 0.5*cos(2*pi*(fc)*n*Ts + phi);  % Rx Signal with a large DC offset

%% Simulation of Costas Adaptive tracking of Carrier Phase
% Moving avergae is computationally faster option than LPF
%
% Assuming the algorithm is running after accumulating 100 samples of $r[n] = s[n] + w[n]$
%
% Initialize $\theta[k]$ to a random number $\in [-\pi/2, \pi/2]$
%
% update step size: $\alpha = 0.005$
%
% $\theta[k]$ gets updated according to following equation:
%
% $\theta[k+1] = \theta[k] - \alpha.avg \{ LPF \{ r[n].cos_{lo}[n] \}.LPF \{ r[n].sin_{lo}[n] \} \} $
avg_len = 1000;
for n_runs=1:100
%% Add AWGN Noise to the received signal so that the $SNR = -10dB$
% $r[n] = s[n] + w[n]$
%
% generate Noise with different power for each simulation
    SNR = 10*rand() - 10;
    r_n = awgn(s_n, SNR);
    % Initialize $\theta[k]$ to be some random number between $[-\pi/2, \pi/2]$
    theta=zeros(1,length(n)/avg_len + 1);
    theta(1) = pi*rand()-0.5*pi;
    % update step size: $\alpha = 0.005$
    alpha_step = 0.075;
    % $\theta[k]$ gets updated according to following equation:
    % $\theta[k+1] = \theta[k] - \alpha.avg{LPF{r[n].cos_{lo}[n]}.LPF{r[n].sin_{lo}}}$
    for k=1:length(n)/avg_len
        cos_lo = cos(2*pi*fc*n((avg_len*k-avg_len+1):avg_len*k)*Ts + theta(k));
        sin_lo = sin(2*pi*fc*n((avg_len*k-avg_len+1):avg_len*k)*Ts + theta(k));
        %J_SD(k+1) = mean(lowpass((r_n(avg_len*k-avg_len+1:avg_len*k).*cos_lo), 0.1).^2);
        theta(k+1) = theta(k) - alpha_step*(mean(r_n(avg_len*k-avg_len+1:avg_len*k).*cos_lo).*mean(r_n(avg_len*k-avg_len+1:avg_len*k).*sin_lo));
    end
    fprintf("Simulation Run %d: SNR: %f theta Init: %f Estimate of Carrier Phase: %f\n",n_runs, SNR, theta(1), theta(end));
    figure(1);set(gcf,'Position',[0,0,2500,1000]);
    plot(theta);hold on;
end
line([0, length(theta)], [-pi/2+phi, -pi/2+phi], 'Color', [1,1,0]);
line([0, length(theta)], [pi/2+phi, pi/2+phi], 'Color', [1,1,0]);
line([0, length(theta)], [phi, phi], 'Color', [1,0,0]);
hold off;

