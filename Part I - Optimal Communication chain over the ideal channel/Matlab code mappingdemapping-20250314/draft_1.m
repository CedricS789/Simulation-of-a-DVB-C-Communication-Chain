%% clc; clear; close all;

%% ----------------- Parameters -----------------
Nbt = 100;                  % Total number of bits
Nbs = 4;                    % Number of bits per symbol
bits_tx = randi([0 1], Nbt, 1); % Random bit sequence
M = 4;                      % Oversampling factor (samples per symbol)
sym_rate = 5e6;             % Symbol rate [5 Msymb/s]
Tsym = 1 / sym_rate;        % Symbol duration
Fs = M * 2*sym_rate;        % Sampling frequency
BWD_RCC = 6e6;              % Bandwidth of the RRC filter
rolloff = 0.3;              % Roll-off factor (β)
span = 6;                   % Filter span in symbols
taps = span * M + 1;        % Total number of taps in the filter

% Define desired frequency shifts (in Hz)
freq_shifts = [-15e6, 0, 15e6]; % Example: -15 MHz, 0 MHz (baseband), and 15 MHz

%% ----------------- Step 1: Random bit generation and bit mapping -----------------
if Nbs > 1
    symb_tx = mapping(bits_tx, Nbs, 'qam');
else
    symb_tx = mapping(bits_tx, Nbs, 'pam');
end

%% ----------------- Plot: QAM Constellation Diagram -----------------
figure;
plot(real(symb_tx), imag(symb_tx), 'o');
title('QAM Constellation Diagram');
xlabel('In-phase');
ylabel('Quadrature');
grid on;
axis equal;

%% ----------------- Step 2: Generate RRC Filter -----------------
% Generate time vector for filter taps
t = (-span/2 : 1/M : span/2) * Tsym;

% Generate RRC filter coefficients
rrc_filter = rcosdesign(rolloff, span, M, 'sqrt');

% Plot the filter impulse response
figure;
plot(t, rrc_filter);
title('Root Raised Cosine Filter Impulse Response');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Plot frequency response
[H, f] = freqz(rrc_filter, 1, 1024, Fs);
figure;
plot(f/1e6, 20*log10(abs(H)));
title('Root Raised Cosine Filter Frequency Response');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
grid on;
xlim([0 Fs/2e6]);

%% ----------------- Step 3: Upsample and filter the signal -----------------
% Upsample the symbols
symbols_up = upsample(symb_tx, M);

% Apply RRC filter to the signal
tx_filtered = filter(rrc_filter, 1, symbols_up);

% Plot time domain signal
t_sig = (0:length(tx_filtered)-1) / Fs;
figure;
subplot(2,1,1);
plot(t_sig*1e6, real(tx_filtered));
title('Real Part of Filtered Signal (Baseband)');
xlabel('Time (μs)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(t_sig*1e6, imag(tx_filtered));
title('Imaginary Part of Filtered Signal (Baseband)');
xlabel('Time (μs)');
ylabel('Amplitude');
grid on;

% Plot spectrum of the baseband filtered signal
N_fft = 2^nextpow2(length(tx_filtered));
TX_freq = fftshift(fft(tx_filtered, N_fft))/length(tx_filtered);
f_axis = (-N_fft/2:N_fft/2-1)*Fs/N_fft;

figure;
plot(f_axis/1e6, 20*log10(abs(TX_freq)));
title('Spectrum of Filtered Signal (Baseband)');
xlabel('Frequency (MHz)');
ylabel('Power Spectral Density (dB)');
grid on;
xlim([-Fs/2e6 Fs/2e6]);

%% ----------------- Step 4: Frequency Shifting -----------------
% Create a figure for the combined spectrum
figure;
hold on;

% Process each frequency shift
for i = 1:length(freq_shifts)
    freq_shift = freq_shifts(i);
    
    % Generate time vector for the entire signal duration
    t_full = (0:length(tx_filtered)-1)' / Fs;
    
    % Apply frequency shift using complex exponential
    tx_shifted = tx_filtered .* exp(1j*2*pi*freq_shift*t_full);
    
    % Plot time domain of frequency-shifted signal (optional)
    if freq_shift ~= 0 % Only plot for non-baseband shifts
        figure;
        subplot(2,1,1);
        plot(t_sig*1e6, real(tx_shifted));
        title(['Real Part of Signal Shifted to ' num2str(freq_shift/1e6) ' MHz']);
        xlabel('Time (μs)');
        ylabel('Amplitude');
        grid on;
        
        subplot(2,1,2);
        plot(t_sig*1e6, imag(tx_shifted));
        title(['Imaginary Part of Signal Shifted to ' num2str(freq_shift/1e6) ' MHz']);
        xlabel('Time (μs)');
        ylabel('Amplitude');
        grid on;
    end
    
    % Compute and plot spectrum of the frequency-shifted signal
    TX_shifted_freq = fftshift(fft(tx_shifted, N_fft))/length(tx_shifted);
    
    % Add to the combined spectrum plot
    plot(f_axis/1e6, 20*log10(abs(TX_shifted_freq)));
end

% Finalize the combined spectrum plot
title('Spectrum of Signals at Different Frequency Shifts');
xlabel('Frequency (MHz)');
ylabel('Power Spectral Density (dB)');
grid on;
xlim([-Fs/2e6 Fs/2e6]);
legend(['Baseband', arrayfun(@(x) [num2str(x/1e6) ' MHz'], freq_shifts(freq_shifts~=0), 'UniformOutput', false)]);
hold off;

% Plot individual shifted spectrums on separate figures
for i = 1:length(freq_shifts)
    if freq_shifts(i) == 0
        continue; % Skip baseband (already plotted)
    end
    
    freq_shift = freq_shifts(i);
    t_full = (0:length(tx_filtered)-1)' / Fs;
    tx_shifted = tx_filtered .* exp(1j*2*pi*freq_shift*t_full);
    
    TX_shifted_freq = fftshift(fft(tx_shifted, N_fft))/length(tx_shifted);
    
    figure;
    plot(f_axis/1e6, 20*log10(abs(TX_shifted_freq)));
    title(['Spectrum of Signal Shifted to ' num2str(freq_shift/1e6) ' MHz']);
    xlabel('Frequency (MHz)');
    ylabel('Power Spectral Density (dB)');
    grid on;
    xlim([-Fs/2e6 Fs/2e6]);
end