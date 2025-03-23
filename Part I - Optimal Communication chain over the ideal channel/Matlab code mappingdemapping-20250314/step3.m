clear; close all; clc;

%% Parameters

Nbps = 2;           % Bits per symbol (256-QAM)
M = Nbps*2^18;           % Number of bits
F = 5e6;            % Symbol frequency
Tbit = 1/(F*Nbps);  % Duration of 1 bit
T = Tbit * M;       % Total duration
OSF = 4;            % Oversampling factor
beta = 0.2;         % Rolloff factor
BW = 6e6;           % Bandwidth
NbTaps = 101;       % Filter taps
Fs = BW*OSF;        % Sampling frequency

%% BER computation parameters

EbN0max = 30;         % in dB
EbN0_domain = 0:1:EbN0max;
iteration = 1;        % averaging BER
BER_array = zeros(1,iteration);
BER = zeros(1,length(EbN0_domain));

%% Signal processing

% Generate random bit sequence
bit_tx = randi([0,1], 1, M);

% QAM Mapping & Upsampling
symb_tx = mapping(bit_tx, Nbps, 'qam');
symb_tx_up = UpS(symb_tx, OSF).';

% Halt Root Raised Cosine (RRC) Filtering at transmitter
filter = RRC_filter(beta, BW, OSF, NbTaps);
signal_tx = conv(symb_tx_up, filter, 'same');

% Signal Power Calculation
signalPowerBaseBand = mean(abs(signal_tx).^2);  % Baseband signal power
signalPowerWideBand = signalPowerBaseBand/2;    % Wideband signal power
Eb = signalPowerWideBand*Tbit;                  % Energy per bit

EbN0_index = 1;
for EbN0 = EbN0_domain
    for i=1:iteration   

        N0 = Eb/(10^(EbN0/10));             % Noise PSD (Power Spectral Density)
        noisePower = N0*F*OSF;              % AWGN power

        % Generate AWGN Noise
        n = sqrt(noisePower)*(randn(1,length(signal_tx)) + 1i*randn(1, length(signal_tx)));

        % Add Noise to Signal
        sn = signal_tx + n.';

        % Half Root Raised Cosine (RRC) Filtering at receptor
        signal_rx = conv(sn, filter, 'same');

        % Downsampling & Demapping
        symb_rx = DownS(signal_rx, OSF).';
        bit_rx = demapping(symb_rx, Nbps, 'qam')';

        % Computing BER
        BER_array(i) = sum(bit_rx ~= bit_tx)/M;
    end
    BER(EbN0_index) = sum(BER_array(:))/iteration;
    EbN0_index = EbN0_index + 1;
end

% Time Vector for Visualization
t = linspace(0, T, M * OSF);
tx_msg = repelem(bit_tx, OSF);
rx_msg = repelem(bit_rx, OSF);

%% Plots

% Dock Figures
set(0, 'DefaultFigureWindowStyle', 'docked');

figure();
semilogy(EbN0_domain,BER(:));hold on;
semilogy(EbN0_domain, berawgn(EbN0_domain,'qam',2^Nbps))
xlabel('Eb/N0 [dB]');
ylabel('BER');
title('Bit Error Rate');
legend('Experimental', 'Theoritical');
ylim([10^-5 1])