clc;
clear all;
close all;

%% =======================================================================
%  DIGITAL FILTER DESIGN USING Z-TRANSFORM IMPLEMENTATION
%  CHEBYSHEV (LOW-PASS) & HIGH-PASS FILTERS (WITHOUT BUILT-IN FUNCTIONS)
%  ENHANCED VERSION WITH VISUAL & ANALYTICAL IMPROVEMENTS
% =======================================================================

Fs = 1000;                  % Sampling frequency (Hz)
t = 0:1/Fs:1;               % Time vector (1 second)
N = length(t);

% Composite test signal (mixture of low and high frequency)
x = sin(2*pi*50*t) + 0.5*sin(2*pi*250*t) + 0.2*randn(size(t));

% SNR levels for analysis
SNR_values = [30 20 10];  

% ------------------------------------------------------------------------
%  CUSTOM FUNCTIONS
% ------------------------------------------------------------------------

% 1️⃣ Manual Gaussian noise adder (no built-in awgn)
manual_awgn = @(sig,SNR_dB) sig + sqrt(var(sig)/(10^(SNR_dB/10))) * randn(size(sig));

% 2️⃣ Compute manual SNR (in dB)
compute_snr = @(orig,noisy) 10*log10(sum(orig.^2)/sum((orig-noisy).^2));

% 3️⃣ Normalize function for plotting
normalize = @(s) s/max(abs(s));

% ------------------------------------------------------------------------
%  EXPERIMENT 1: CHEBYSHEV FILTER (LOW-PASS)
% ------------------------------------------------------------------------
b0 = 0.2066; b1 = 0.4131; b2 = 0.2066;
a1 = -0.3695; a2 = 0.1958;

figure('Name','Experiment 1 - Chebyshev Filter (Z-transform)','NumberTitle','off');
colormap turbo;

for i = 1:length(SNR_values)
    % Step 1: Add Gaussian noise manually
    noisy = manual_awgn(x, SNR_values(i));

    % Step 2: Apply Chebyshev filter recursively
    y = zeros(size(noisy));
    for n = 3:length(noisy)
        y(n) = b0*noisy(n) + b1*noisy(n-1) + b2*noisy(n-2) ...
             - a1*y(n-1) - a2*y(n-2);
    end

    % Step 3: Compute SNR improvement
    input_snr = compute_snr(x, noisy);
    output_snr = compute_snr(x, y);

    % Step 4: Time Domain
    subplot(length(SNR_values),2,2*i-1);
    plot(t,normalize(noisy),'r','LineWidth',1);
    hold on; plot(t,normalize(y),'b','LineWidth',1.2);
    title(sprintf('Time Domain @ %d dB | SNR_{out}=%.2f dB', SNR_values(i), output_snr));
    xlabel('Time (s)'); ylabel('Normalized Amplitude');
    legend('Noisy','Filtered'); grid on;

    % Step 5: Frequency Spectrum
    subplot(length(SNR_values),2,2*i);
    f = (0:N-1)*(Fs/N);
    Y = abs(fft(y));
    plot(f,normalize(Y),'b','LineWidth',1.3);
    title(['Frequency Spectrum (Chebyshev LPF) @ ' num2str(SNR_values(i)) ' dB']);
    xlabel('Frequency (Hz)'); ylabel('|H(f)|'); grid on;
end

% ------------------------------------------------------------------------
%  EXPERIMENT 2: HIGH-PASS FILTER (Z-transform)
% ------------------------------------------------------------------------
a = 0.9; % smoothing factor (closer to 1 → sharper HPF)

figure('Name','Experiment 2 - High-Pass Filter (Z-transform)','NumberTitle','off');
for i = 1:length(SNR_values)
    noisy = manual_awgn(x, SNR_values(i));
    y_hp = zeros(size(noisy));
    for n = 2:length(noisy)
        y_hp(n) = a*(y_hp(n-1) + noisy(n) - noisy(n-1));
    end
    
    % Compute output SNR
    output_snr = compute_snr(x, y_hp);

    % Time domain
    subplot(length(SNR_values),2,2*i-1);
    plot(t,normalize(noisy),'r','LineWidth',1);
    hold on; plot(t,normalize(y_hp),'b','LineWidth',1.2);
    title(sprintf('High-Pass Output @ %d dB | SNR_{out}=%.2f dB', SNR_values(i), output_snr));
    xlabel('Time (s)'); ylabel('Normalized Amplitude');
    legend('Noisy','Filtered'); grid on;

    % Frequency domain
    subplot(length(SNR_values),2,2*i);
    f = (0:N-1)*(Fs/N);
    Yhp = abs(fft(y_hp));
    plot(f,normalize(Yhp),'b','LineWidth',1.3);
    title(['Frequency Spectrum (High-Pass) @ ' num2str(SNR_values(i)) ' dB']);
    xlabel('Frequency (Hz)'); ylabel('|H(f)|'); grid on;
end

% ------------------------------------------------------------------------
%  FREQUENCY RESPONSE (Z-Domain Theoretical)
% ------------------------------------------------------------------------
w = linspace(0,pi,512);
H_cheby = (b0 + b1*exp(-1j*w) + b2*exp(-2j*w)) ./ (1 + a1*exp(-1j*w) + a2*exp(-2j*w));
H_hp = a * (1 - exp(-1j*w)) ./ (1 - a*exp(-1j*w));

figure('Name','Z-Transform Frequency Responses','NumberTitle','off');
subplot(2,1,1);
plot(w/pi, abs(H_cheby),'LineWidth',1.5);
title('Chebyshev Low-Pass |H(e^{j\omega})|');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude'); grid on;

subplot(2,1,2);
plot(w/pi, abs(H_hp),'LineWidth',1.5);
title('High-Pass |H(e^{j\omega})|');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude'); grid on;

% ------------------------------------------------------------------------
%  ENERGY ANALYSIS: Before vs After Filtering
% ------------------------------------------------------------------------
figure('Name','Energy Distribution Comparison','NumberTitle','off');
E_in = sum(x.^2);
E_cheby = sum(y.^2);
E_hp = sum(y_hp.^2);
bar([E_in E_cheby E_hp]);
set(gca,'XTickLabel',{'Original','Chebyshev','High-Pass'});
ylabel('Signal Energy'); title('Energy Comparison After Filtering');
grid on;

% ------------------------------------------------------------------------
%  SUMMARY DISPLAY
% ------------------------------------------------------------------------
disp('===============================================================');
disp(' ENHANCED DIGITAL FILTER DESIGN USING Z-TRANSFORM ');
disp(' Implemented Filters: Chebyshev (Low-pass) & High-pass ');
disp(' Added Features: Dynamic SNR, Spectrum, and Energy Analysis ');
disp(' No built-in filtering or AWGN functions used ');
disp('===============================================================');
