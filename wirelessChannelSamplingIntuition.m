% Frequency Selective Fading Channel - Coh BW << Sys BW
sys_bw = 5e6;         % 5MHz System Bandwidth
fs = 20e6;            % 20MHz Sampling Frequency
coh_bw = 500e3;       % 500KHz Coherence Bandwidth
delay_spread = 1/2/coh_bw;
numPaths = 21;

path_delays = rand(1,numPaths)*delay_spread;
%path_delays = [4.29006834809595e-07,8.00334711846513e-07,4.68572355657067e-07,4.08755089279763e-08,5.88215931860844e-07,7.38653490429940e-08,5.91447482737160e-07,5.00922207668102e-07,2.54891877514721e-07,4.75009968485334e-07,2.10720117613309e-07,6.08287445617845e-07,9.15735132287457e-07,5.47016878644843e-07,4.66139552549696e-07,1.54551725249510e-07,6.40689499329106e-07,2.71760535229616e-07,8.62821600770533e-07,3.66709073895298e-07,4.21534639554125e-07];
path_gains = raylrnd(1, 1, numPaths);

t = (-200/fs):1/fs:(200/fs);
channel_analog = zeros(size(t));
for i=1:numPaths
    kernel(i, :) = path_gains(i)*sinc(fs*(t-path_delays(i)));
    channel_analog = channel_analog + kernel(i,:);
end
sample_points = 0:1/fs:2*delay_spread;
channel_taps = interp1(t, channel_analog, sample_points, 'spline');

%% For plotting lets sample the sinc at more points
t = (-200/fs):1/10/fs:(200/fs);
channel_analog_up = zeros(size(t));
figure(1);
stem(path_delays*1e9, path_gains, 'g--', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');hold on;
for i=1:numPaths
    kernel_up(i, :) = path_gains(i)*sinc(fs*(t-path_delays(i)));
    plot((t)*1e9, kernel_up(i, :), ':', 'lineWidth', 1);
    channel_analog_up = channel_analog_up + kernel_up(i,:);
end
plot(t*1e9, channel_analog_up, 'b.-', 'lineWidth', 1);
xline(0, 'k')
xlabel('nsec');
ylabel('RF Channel(time-domain)');
xlim([-1/fs*1e9, 2*max(path_delays)*1e9]);
title('Channel Visualization in Analog Time-Domain');


figure(2);
% pwelch(channel_analog, [], [], [], fs);hold on;
% pwelch(kernel1, [], [], [], fs);hold off;
channel_freq = fftshift(fft(channel_analog));
system_freq = fftshift(fft(sinc(fs*(-200/fs:1/fs:200/fs))));
plot(linspace(-fs/2, fs/2, numel(channel_analog))/1e6, 20*log10(abs(channel_freq)), 'b', 'lineWidth', 1.5);hold on;
plot(linspace(-fs/2, fs/2, numel(channel_analog))/1e6, 20*log10(abs(system_freq)), 'r', 'lineWidth', 1.5);
title('Channel Visualization in Frequency Domain');
xlabel('Freq(MHz)');
ylabel('Channel Magnitude Response(dB)');

figure(3);
stem(channel_taps, '^', 'lineWidth', 1, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
xlabel('Taps');
ylabel('Channel Taps, h_l[n]');