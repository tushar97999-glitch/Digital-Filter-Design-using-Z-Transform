clc;
clear all;
close all;

%% =======================================================================
%  DIGITAL FILTER DESIGN USING Z-TRANSFORM IMPLEMENTATION
%  FOR CHEBYSHEV AND HIGH-PASS FILTERS (WITHOUT BUILT-IN FUNCTIONS)
%  ENHANCED VERSION WITH EXTRA OUTPUTS & STATISTICS
% =======================================================================

Fs = 1000;                % Sampling frequency (Hz)
t = 0:1/Fs:1;             % Time vector (1 second)
x = randn(size(t));       % Random signal

SNR_values = [30 20 10];  % Gaussian noise levels (in dB)

% -----------------------------------------------------------------------
%  FUNCTION DEFINITIONS
% -----------------------------------------------------------------------

% (1) Manual function to add Gaussian noise
manual_awgn = @(sig,SNR_dB) ...
    sig + sqrt(var(sig)/(10^(SNR_dB/10))) * randn(size(sig));

% (2) Manual SNR computation
manual_snr = @(orig,noisy) 10*log10(sum(orig.^2)/sum((orig-noisy).^2));

% (3) Manual RMS computation
manual_rms = @(sig) sqrt(mean(sig.^2));

%% =======================================================================
%  EXPERIMENT 1: CHEBYSHEV FILTER (LOW-PASS)
% =======================================================================
b0 = 0.2066;
b1 = 0.4131;
b2 = 0.2066;
a1 = -0.3695;
a2 = 0.1958;

figure('Name','Experiment 1 - Chebyshev Filter (Z-transform)');

disp('---------------------------------------------------------------');
disp(' CHEBYSHEV FILTER PERFORMANCE (Manual Implementation)');
disp('---------------------------------------------------------------');
fprintf('   SNR_in(dB) | SNR_out(dB) | Energy_in | Energy_out | Mean_out | RMS_out\n');
disp('--------------------------------------------------------------------------');

for i = 1:length(SNR_values)
    % Step 1: Add noise
    noisy = manual_awgn(x, SNR_values(i));
    
    % Step 2: Manual recursive filtering
    y = zeros(size(noisy));
    for n = 3:length(noisy)
        y(n) = b0*noisy(n) + b1*noisy(n-1) + b2*noisy(n-2) ...
             - a1*y(n-1) - a2*y(n-2);
    end

    % Step 3: Signal metrics
    SNR_in = manual_snr(x, noisy);
    SNR_out = manual_snr(x, y);
    Energy_in = sum(noisy.^2);
    Energy_out = sum(y.^2);
    Mean_out = mean(y);
    RMS_out = manual_rms(y);

    % Print to console
    fprintf('      %4.1f     |    %4.2f     |  %7.2f  |  %7.2f  |  %7.3f  |  %7.3f\n',...
        SNR_in, SNR_out, Energy_in, Energy_out, Mean_out, RMS_out);

    % Step 4: Plot in time domain
    subplot(length(SNR_values),3,3*i-2);
    plot(t,noisy,'r'); title(['Noisy (' num2str(SNR_values(i)) ' dB)']);
    xlabel('Time (s)'); ylabel('Amplitude'); grid on;

    subplot(length(SNR_values),3,3*i-1);
    plot(t,y,'b'); title('Filtered Output');
    xlabel('Time (s)'); ylabel('Amplitude'); grid on;

    % Step 5: Manual Power Spectral Density (FFT)
    N = length(noisy);
    f = (0:N-1)*(Fs/N);
    PSD_in = abs(fft(noisy)).^2/N;
    PSD_out = abs(fft(y)).^2/N;

    subplot(length(SNR_values),3,3*i);
    plot(f,10*log10(PSD_in/max(PSD_in)),'r'); hold on;
    plot(f,10*log10(PSD_out/max(PSD_out)),'b');
    title('Normalized PSD (dB)');
    xlabel('Frequency (Hz)'); ylabel('Power (dB)');
    legend('Noisy','Filtered'); grid on;
end

%% =======================================================================
%  EXPERIMENT 2: HIGH-PASS FILTER (Z-transform)
% =======================================================================
a = 0.9; % controls cutoff frequency

figure('Name','Experiment 2 - High-Pass Filter (Z-transform)');

disp(' ');
disp('---------------------------------------------------------------');
disp(' HIGH-PASS FILTER PERFORMANCE (Manual Implementation)');
disp('---------------------------------------------------------------');
fprintf('   SNR_in(dB) | SNR_out(dB) | Energy_in | Energy_out | Mean_out | RMS_out\n');
disp('--------------------------------------------------------------------------');

for i = 1:length(SNR_values)
    % Step 1: Add noise
    noisy = manual_awgn(x, SNR_values(i));
    
    % Step 2: Apply HPF using difference equation
    y_hp = zeros(size(noisy));
    for n = 2:length(noisy)
        y_hp(n) = a*(y_hp(n-1) + noisy(n) - noisy(n-1));
    end

    % Step 3: Signal metrics
    SNR_in = manual_snr(x, noisy);
    SNR_out = manual_snr(x, y_hp);
    Energy_in = sum(noisy.^2);
    Energy_out = sum(y_hp.^2);
    Mean_out = mean(y_hp);
    RMS_out = manual_rms(y_hp);

    fprintf('      %4.1f     |    %4.2f     |  %7.2f  |  %7.2f  |  %7.3f  |  %7.3f\n',...
        SNR_in, SNR_out, Energy_in, Energy_out, Mean_out, RMS_out);

    % Step 4a: Time Domain
    subplot(length(SNR_values),3,3*i-2);
    plot(t,noisy,'r'); title(['Noisy (' num2str(SNR_values(i)) ' dB)']);
    xlabel('Time (s)'); ylabel('Amplitude'); grid on;

    subplot(length(SNR_values),3,3*i-1);
    plot(t,y_hp,'b'); title('Filtered Output');
    xlabel('Time (s)'); ylabel('Amplitude'); grid on;

    % Step 4b: Frequency Domain (Manual FFT)
    Nfft = length(noisy);
    f = (0:Nfft-1)*(Fs/Nfft);
    Y = abs(fft(y_hp));
    subplot(length(SNR_values),3,3*i);
    plot(f,Y,'b'); title('FFT Magnitude');
    xlabel('Frequency (Hz)'); ylabel('|Y(f)|'); grid on;
end

%% =======================================================================
%  THEORETICAL FREQUENCY RESPONSES (Z-DOMAIN)
% =======================================================================
w = linspace(0,pi,512);

H_cheby = (b0 + b1*exp(-1j*w) + b2*exp(-2j*w)) ./ ...
          (1 + a1*exp(-1j*w) + a2*exp(-2j*w));

H_hp = a * (1 - exp(-1j*w)) ./ (1 - a*exp(-1j*w));

figure('Name','Theoretical Frequency Responses');
subplot(2,1,1);
plot(w/pi, abs(H_cheby),'LineWidth',1.5);
title('Chebyshev Filter |H(e^{j\omega})|');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude'); grid on;

subplot(2,1,2);
plot(w/pi, abs(H_hp),'LineWidth',1.5);
title('High-Pass Filter |H(e^{j\omega})|');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude'); grid on;

%% =======================================================================
%  DISPLAY FINAL SUMMARY
% =======================================================================
disp(' ');
disp('===============================================================');
disp(' DIGITAL FILTER DESIGN USING Z-TRANSFORM IMPLEMENTATION ');
disp(' Implemented Filters: Chebyshev (Low-pass) and High-pass ');
disp(' Added Outputs: SNR, Energy, PSD, Mean, RMS ');
disp(' No built-in filter or awgn functions used ');
disp('===============================================================');
