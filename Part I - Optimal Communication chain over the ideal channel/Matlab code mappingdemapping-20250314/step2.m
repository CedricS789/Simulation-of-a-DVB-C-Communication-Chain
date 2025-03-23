clear; close all; clc;

%% Parameters

Nbps = 4;           % Bits per symbol (256-QAM)
M = Nbps*100;       % Number of bits
F = 5e6;            % Symbol frequency
OSF = 4;            % Oversampling factor
%Fs = F * OSF;      % Sampling frequency
T = M / (F * Nbps); % Total duration
beta = 0.2;         % Rolloff factor
BW = 6e6;           % Bandwidth
NbTaps = 101;       % Filter taps
Fs = BW * OSF;      % Sampling frequency

%% Signal processing

% Generate random bit sequence
bit_tx = randi([0,1], 1, M);

% QAM Mapping & Upsampling
symb_tx = mapping(bit_tx, Nbps, 'qam');
symb_tx_up = UpS(symb_tx, OSF).';

% Root Raised Cosine (RRC) Filtering
filter = RRC_filter(beta, BW, OSF, NbTaps);
%signal_tx = conv(symb_tx_up, filter, 'same');
%signal_rx = conv(signal_tx, filter, 'same');
signal_tx = conv(symb_tx_up, filter, 'full');                               % Full convolution
signal_tx = signal_tx(floor(length(filter)/2) + (1:length(symb_tx_up)));    % Remove delay

signal_rx = conv(signal_tx, filter, 'full');                                % Full convolution
signal_rx = signal_rx(floor(length(filter)/2) + (1:length(signal_tx)));     % Remove delay

% Downsampling & Demapping
symb_rx = DownS(signal_rx, OSF).';
bit_rx = demapping(symb_rx, Nbps, 'qam')';

% Time Vector for Visualization
t = linspace(0, T, M * OSF);
tx_msg = repelem(bit_tx, OSF);
rx_msg = repelem(bit_rx, OSF);

% Compute Bit Error
error = norm(bit_rx - bit_tx) / M;

%% Plots

% Dock Figures
set(0, 'DefaultFigureWindowStyle', 'docked');

% Transmitted Symbols Visualization
figure;
subplot(2,1,1);
scatter(real(symb_tx), imag(symb_tx), 20, 'o');
hold on;
plot([0 0], [-1.5 1.5], 'k', 'LineWidth', 1);
plot([-1.5 1.5], [0 0], 'k', 'LineWidth', 1);
grid on;
axis equal;
xlabel('Re'); ylabel('Im');
title('Transmitted Symbols');

subplot(2,1,2);
plot(t,tx_msg);
xlabel('Time'); ylabel('Bit Value');
title('Transmitted Bit Stream');
axis([t(1) t(end) -0.5 1.5])

% Received Symbols Visualization
figure;
subplot(2,1,1);
scatter(real(symb_rx), imag(symb_rx), 20, 'o');
hold on;
plot([0 0], [-1.5 1.5], 'k', 'LineWidth', 1);
plot([-1.5 1.5], [0 0], 'k', 'LineWidth', 1);
grid on;
axis equal;
xlabel('Re'); ylabel('Im');
title('Received Symbols');

subplot(2,1,2);
plot(t,rx_msg);
xlabel('Time'); ylabel('Bit Value');
title('Received Bit Stream');
axis([t(1) t(end) -0.5 1.5])
