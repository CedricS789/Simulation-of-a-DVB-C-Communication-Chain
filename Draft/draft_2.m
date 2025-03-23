%% Clear workspace and close all figures
clc; clear; close all;

%% ----------------- Parameters -----------------
Nbt = 1000;                              % Total number of bits to transmit
Nbs = 4;                                % Number of bits per symbol (e.g., 4 for 16-QAM)
bits_tx = randi([0 1], Nbt, 1);         % Generate random bit sequence of 0s and 1s

M = 4;                                  % Oversampling factor (samples per symbol)
sym_rate = 5e6;                         % Symbol rate in symbols per second [5 Msymb/s]
Tsym = 1 / sym_rate;                    % Symbol duration in seconds
Fs = M * sym_rate;                      % Sampling frequency in Hz (samples per second)

rolloff = 0.2;                          % Roll-off factor for RRC filter
bandwidth = (1 + rolloff) * sym_rate;   % Bandwidth of the filter, based on roll-off and symbol rate
span = 50;                              % Filter span in number of symbol periods
taps = span * M + 1;                    % Total number of filter taps (including center tap)

%% ----------------- Step 1: Random Bit Generation and Symbol Mapping -----------------
if Nbs > 1
    symb_tx = mapping(bits_tx, Nbs, 'qam');  % QAM mapping for Nbs > 1
else
    symb_tx = mapping(bits_tx, Nbs, 'pam');  % PAM mapping for Nbs = 1
end

%% ----------------- Step 2: Generate and Apply RRC Filter -----------------
symbols_up = upsample(symb_tx, M);                  % Upsampled symbol sequence
rrc_filter = rcosdesign(rolloff, span, M, 'sqrt');  % RRC filter with square-root response
tx_filtered = filter(rrc_filter, 1, symbols_up);    % Filtered transmit signal

%% ----------------- Print Relevant Parameters -----------------
if Nbs > 1
    modulation_type = 'QAM';             % Quadrature Amplitude Modulation
else
    modulation_type = 'PAM';             % Pulse Amplitude Modulation
end

fprintf('Total number of bits: %d\n', Nbt);
fprintf('Bits per symbol: %d\n', Nbs);
fprintf('Modulation type: %s\n', modulation_type);
fprintf('Symbol rate: %.2f Msymbols/s\n', sym_rate / 1e6);
fprintf('Symbol period duration: %.2f μs\n', Tsym * 1e6);
fprintf('Oversampling factor: %d\n', M);
fprintf('Sampling frequency: %.2f MHz\n', Fs / 1e6);
fprintf('Roll-off factor: %.2f\n', rolloff);
fprintf('Filter span: %d symbols\n', span);
fprintf('Number of filter taps: %d\n', taps);
fprintf('Bandwidth: %.2f MHz\n', bandwidth / 1e6);

%% ----------------- Plot: RRC Filter Impulse Response with Symbol Rate Sampling -----------------
t = (-span*M/2 : 1 : span*M/2) / Fs;  % Time in seconds, centered around 0

% Create a full raised cosine filter (not square root) to demonstrate Nyquist property
% Full Nyquist filter = RRC * RRC (convolution in time domain)
rc_filter = rcosdesign(rolloff, span, M, 'normal');  % Normal raised cosine filter

% Symbol-rate sampling points (every M samples)
sample_indices = 1:M:length(rc_filter);
sampled_times = t(sample_indices);
sampled_values = rc_filter(sample_indices);

figure;
plot(t, rc_filter, 'b-', 'LineWidth', 1.5);  % Plot the continuous filter response
hold on;
plot(sampled_times, sampled_values, 'ro', 'MarkerSize', 8, 'LineWidth', 2);  % Plot the sampled points
title('Raised Cosine Filter Impulse Response');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
legend('Continuous Filter Response', 'Symbol-Rate Samples');
hold off;

% Calculate and display the maximum deviation of non-center samples from zero
center_idx = ceil(length(sampled_times)/2);
nonzero_samples = sampled_values;
nonzero_samples(center_idx) = [];  % Remove the center sample
max_deviation = max(abs(nonzero_samples));
fprintf('Maximum deviation of non-center samples from zero: %.8f\n', max_deviation);

%% ----------------- Plot: RRC Filter Frequency Response -----------------
N = 1024;                               % Number of FFT points for better resolution
H = fft(rrc_filter, N);                 % Compute FFT of the filter with zero-padding
H_shifted = fftshift(H);                % Shift FFT to center zero frequency
f = (-N/2 : N/2 - 1) * (Fs/N);          % Frequency axis in Hz
figure;
plot(f / 1e6, 20 * log10(abs(H_shifted)));
title('Root Raised Cosine Filter Frequency Response');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
grid on;

%% ----------------- Plot: QAM Constellation Diagram -----------------
figure;
plot(real(symb_tx), imag(symb_tx), 'o');    % Scatter plot of real vs. imaginary parts
title('QAM Constellation Diagram');
xlabel('In-phase');
ylabel('Quadrature');
grid on;
axis equal;                                 % Ensure equal scaling for accurate representation

%% ----------------- Plot: Time Domain of the Filtered Signal -----------------
t_sig = (0 : length(tx_filtered) - 1) / Fs;     % Time in seconds
figure;
subplot(2, 1, 1);
plot(t_sig * 1e6, real(tx_filtered));           % Real part, time in microseconds
title('Real Part of Filtered Signal');
xlabel('Time (μs)');
ylabel('Amplitude');
grid on;
subplot(2, 1, 2);
plot(t_sig * 1e6, imag(tx_filtered));           % Imaginary part, time in microseconds
title('Imaginary Part of Filtered Signal');
xlabel('Time (μs)');
ylabel('Amplitude');
grid on;

%% ----------------- Plot: Spectrum of the Filtered Signal -----------------
N_fft = 2^nextpow2(length(tx_filtered));                                % FFT size as next power of 2 for efficiency
TX_freq = fftshift(fft(tx_filtered, N_fft)) / length(tx_filtered);      % Normalized FFT
f_axis = (-N_fft / 2 : N_fft / 2 - 1) * (Fs / N_fft);                   % Frequency axis in Hz
figure;
plot(f_axis / 1e6, 20 * log10(abs(TX_freq)));
title('Spectrum of Filtered Signal');
xlabel('Frequency (MHz)');
ylabel('Power Spectral Density (dB)');
grid on;
xlim([-Fs/2e6, Fs/2e6]);

%% ----------------- Demonstration of ISI Cancellation -----------------
% Create a sequence of symbols with a single '1' and rest zeros to demonstrate impulse response
test_symbols = zeros(span+1, 1);
test_symbols(ceil(length(test_symbols)/2)) = 1;  % Single '1' in the middle

% Upsample and filter with the raised cosine (not square root)
test_upsampled = upsample(test_symbols, M);
test_filtered = filter(rc_filter, 1, test_upsampled);

% Time axis for the filtered signal
t_test = (0:length(test_filtered)-1)/Fs;

% Sample at symbol rate (every M samples)
sample_indices_test = M:M:length(test_filtered);
sampled_times_test = t_test(sample_indices_test);
sampled_values_test = test_filtered(sample_indices_test);

figure;
plot(t_test, test_filtered, 'b-', 'LineWidth', 1.5);
hold on;
plot(sampled_times_test, sampled_values_test, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
title('ISI Cancellation Demonstration (Raised Cosine Response to a Single Symbol)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
legend('Continuous Filter Response', 'Symbol-Rate Samples');
hold off;