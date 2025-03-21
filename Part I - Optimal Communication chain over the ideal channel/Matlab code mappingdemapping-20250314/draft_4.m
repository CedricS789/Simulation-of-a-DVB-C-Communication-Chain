clc; clear; close all;

%% ----------------- Parameters -----------------
Nbt = 600;                              % Total number of bits
Nbs = 4;                                % Number of bits per symbol
bits_tx = randi([0 1], Nbt, 1);         % Random bit sequence

M = 4;                                  % Oversampling factor (samples per symbol)
sym_rate = 5e6;                         % Symbol rate [5 Msymb/s]
Tsym = 1 / sym_rate;                    % Symbol duration
Fs = M * sym_rate;                      % Sampling frequency

rolloff = 0.3;                          % Roll-off factor (β)
span    = 6;                            % Filter span in symbols
taps = span * M + 1;                    % Total number of taps in the filter

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

%% ----------------- Step 2: Root Raised Cosine Filter Design -----------------

% Time vector for the filter
t = (-span*M/2 : 1 : span*M/2) / Fs;

% RRC filter impulse response
h_rrc = rcosdesign(rolloff, span, M, 'sqrt');

% Verify the filter length
if length(h_rrc) ~= taps
    error('RRC filter length does not match expected length.');
end

%% ----------------- Step 3: Upsampling and Filtering -----------------

% Upsample the symbols
symb_tx_upsampled = upsample(symb_tx, M);

% Filter the upsampled symbols with full convolution and trim delay
symb_tx_filtered = conv(symb_tx_upsampled, h_rrc, 'full');
delay = (length(h_rrc) - 1)/2; % Calculate group delay
symb_tx_filtered = symb_tx_filtered(delay+1:end-delay); % Trim delay

%% ----------------- Plot: RRC Filter Impulse Response -----------------
figure;
plot(t, h_rrc);
title('Root Raised Cosine Filter Impulse Response');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% ----------------- Plot: Upsampled Signal -----------------
figure;
plot(1:length(symb_tx_upsampled), real(symb_tx_upsampled));
title('Upsampled Signal');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

%% ----------------- Plot: Filtered Signal -----------------
figure;
plot(1:length(symb_tx_filtered), real(symb_tx_filtered));
title('Filtered Signal');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;
