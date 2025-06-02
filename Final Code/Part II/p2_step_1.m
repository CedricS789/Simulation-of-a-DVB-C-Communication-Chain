%% =================== p2_step_1  - CFO and Phase Offsett ========================
clear; close all; clc;
addpath('../Part I/p1_functions');
addpath('p2_functions')


%% =================== Load Simulation Parameters  ===================
Nbps = 4;
params = initParameters(Nbps);
displayParameters(params);
NumBits = params.timing.NumBits;
ModType = params.modulation.ModulationType;
ModOrder= params.modulation.ModulationOrder;
OSF = params.sampling.OversamplingFactor;
SymRate = params.timing.SymbolRate;
BitRate = params.timing.BitRate;
Fs = params.sampling.SamplingFrequency;
Ts = params.sampling.SamplePeriod;
Beta = params.filter.RolloffFactor;
NumTaps = params.filter.NumFilterTaps;


%% =================== Communication Chain ===================
% ---- CFO Parameters ----
Fc = 600e6;                             % Carrier frequency in Hz
ppm = 1;
delta_cfo_ppm = ppm * 1e-6 * Fc;        % Frequency offset in Hz
delta_omega = 2 * pi * delta_cfo_ppm;   % Frequency offset in rad/s
phi_0 = 0;                              % Phase offset in rad

% --- Transmitter  ---
bit_tx = randi([0, 1], 1, NumBits).';
symb_tx = mapping(bit_tx, Nbps, ModType);
symb_tx_up = upSampler(symb_tx, OSF).';
g_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);
signal_tx = applyFilter(symb_tx_up, g_rrc, NumTaps);
signalPower_tx = mean(abs(signal_tx).^2);
Eb = signalPower_tx / BitRate;

% -- Introduce Noise --
EbN0dB = 1000;
signal_tx_noisy = addAWGN(signal_tx, Eb, EbN0dB, OSF, SymRate);

% --- Introduce CFO and phase offset ---
time_vector = (0 : length(signal_tx) - 1).' * Ts;
signal_tx_distorted  = signal_tx_noisy .* exp(1j * (delta_omega * time_vector + phi_0)); 

% --- Receiver Chain ---
signal_rx = applyFilter(signal_tx_distorted, g_rrc, NumTaps);
symb_rx_down = downSampler(signal_rx, OSF);
bit_rx = demapping(symb_rx_down, Nbps, ModType); 
bit_rx = bit_rx(:);  


%% =================== Generate Plots  ===================
bits_to_plot = min(params.timing.NumBits, 100 * Nbps); 
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx_up, symb_rx_down);
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);
plotPSD_Tx_Rx(signal_tx, signal_rx, Fs);
plotBasebandFrequencyResponse(signal_tx, signal_rx, Fs);