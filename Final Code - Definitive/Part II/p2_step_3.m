%% ====================================== p2_step_3 - Gardner  ======================================
close all; clear; clc;
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
Tsymb = params.timing.SymbolPeriod;
BitRate = params.timing.BitRate;
Fs = params.sampling.SamplingFrequency;
Ts = params.sampling.SamplePeriod;
Beta = params.filter.RolloffFactor;
NumTaps = params.filter.NumFilterTaps;

% --- BER Curve Parameters ---
EbN0_min_dB         = params.simulation.EbN0_min_dB;                % Minimum Eb/N0 value in dB
EbN0_max_dB         = params.simulation.EbN0_max_dB;                % Maximum Eb/N0 value in dB
EbN0_step_dB        = params.simulation.EbN0_step_dB;               % Step size for Eb/N0 sweep in dB
iterations          = params.simulation.iterations_per_EbN0;        % Iterations for averaging BER at each Eb/N0 point
EbN0_domain_dB      = params.simulation.EbN0_domain_dB;             % Range of Eb/N0 values to simulate (dB)
num_EbN0_points     = length(EbN0_domain_dB);  

% ---- CFO and Sample Time Offset Parameters ----
Fc = 600e6;                                                       % Carrier frequency in Hz
ppm = 0;
delta_cfo = ppm * 1e-6 * Fc;                                      % Frequency offset in Hz (0.08 ppm)
phi_0 = 0;                                                        % Phase offset in rad
timing_offset_percent = 0.2;                                      % Normalized timing offset (% of symbol period)
initial_offset_samples = round(timing_offset_percent * OSF);      % Initial offset in samples (1 sample)
kappa = 0.1;


%% =================== Communication Chain ===================
% --- Transmitter  ---
bit_tx = randi([0, 1], 1, NumBits).';
symb_tx = mapping(bit_tx, Nbps, ModType);
symb_tx_up = upSampler(symb_tx, OSF).';
g_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);
signal_tx_filtered = applyFilter(symb_tx_up, g_rrc, NumTaps);
signalPower_tx = mean(abs(signal_tx_filtered).^2);
Eb = signalPower_tx / BitRate;

time_vector = (0 : length(signal_tx_filtered) - 1).' * Ts;

EbN0dB = 1e10;

signal_tx_noisy = addAWGN(signal_tx_filtered, Eb, EbN0dB, OSF, SymRate);

signal_tx_distorted = signal_tx_noisy .* exp(1j * ((2 * pi * delta_cfo) * time_vector + phi_0));

signal_rx_matched_filtered  = applyFilter(signal_tx_distorted, g_rrc, NumTaps);

signal_rx_distorted = circshift(signal_rx_matched_filtered, initial_offset_samples);

symb_rx_down_partial = downSampler(signal_rx_distorted, OSF/2);

[symb_rx_corected_down, time_shift_errors] = gardner(symb_rx_down_partial, kappa, OSF/(OSF/2));

time_vector_symb = (0 : length(symb_rx_corected_down) - 1).' * Tsymb;

symb_rx_down_compensated = symb_rx_corected_down .* exp(-1j * (2 * pi * delta_cfo * time_vector_symb));

bit_rx = demapping(symb_rx_down_compensated, Nbps, ModType); 


%% =================== Generate Plots  ===================
plotConstellation_Tx_Rx(ModOrder, ModType, symb_rx_down_partial, symb_rx_corected_down);
bits_to_plot = min(params.timing.NumBits, 100 * Nbps); 
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);
plotPSD_Tx_Rx(signal_tx_filtered, signal_rx_matched_filtered, Fs);
plotBasebandFrequencyResponse(signal_tx_filtered, signal_rx_matched_filtered, Fs);
plotFilterCharacteristics(g_rrc, Beta, Fs, OSF);