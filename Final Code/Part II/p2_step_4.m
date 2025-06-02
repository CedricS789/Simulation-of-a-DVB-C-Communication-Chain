%% =================== Step 4 - Frame and Frequency Acquisition and Entire Communication Chain ===================
clear; close all; clc;
addpath('../Part I/p1_functions');
addpath('p2_functions');

%% ========================================== Load Simulation Parameters  ==========================================
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

% ---- CFO and sample time offset Parameters ----
Fc = 600e6;
ppm = 2;
delta_cfo = ppm * 1e-6 * Fc;
phi_0  = 0;
timing_offset_percent = 0.1;
initial_offset_samples = round(timing_offset_percent * OSF);
kappa = 0.1;

pilot_position = 50;
averaging_window = 16; 

%% ========================================== Communication Chain ==========================================
% --- Transmitter  ---
bit_tx = randi([0, 1], 1, NumBits).';
symb_tx = mapping(bit_tx, Nbps, ModType);
pilot_size = floor(length(symb_tx) * 0.4);

unuseful = symb_tx(1 : pilot_position - 1);
pilot = symb_tx(pilot_position : pilot_position + pilot_size - 1);
symb_tx_up = upSampler(symb_tx, OSF).';
g_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);
signal_tx_filtered  = applyFilter(symb_tx_up, g_rrc, NumTaps);
signalPower_tx  = mean(abs(signal_tx_filtered).^2);
Eb = signalPower_tx / BitRate;

time_vector = (0 : length(signal_tx_filtered) - 1).' * Ts;
time_vector_symb = (0 : length(symb_tx) - 1).' * Tsymb;

% -- Introduce Noise --
EbN0dB = 1e10;
signal_tx_noisy = addAWGN(signal_tx_filtered, Eb, EbN0dB, OSF, SymRate);

% --- Introduce CFO, phase offset and sample time offset ---
signal_tx_distorted  = signal_tx_noisy .* exp(1j * (2 * pi * delta_cfo * time_vector + phi_0));
signal_tx_distorted = circshift(signal_tx_distorted, initial_offset_samples);

% --- Receiver Chain ---
signal_rx_matched_filtered  = applyFilter(signal_tx_distorted, g_rrc, NumTaps);
symb_rx_down = downSampler(signal_rx_matched_filtered, OSF);
[symb_rx_corected_down, time_shift_errors] = gardner(signal_rx_matched_filtered, kappa, OSF);
[toa, delta_cfo_hat] = frameFreqAcquisition(pilot, symb_rx_corected_down, averaging_window, Tsymb);
symb_rx_corected_down = symb_rx_corected_down .* exp(-1j * (2 * pi * delta_cfo_hat * time_vector_symb));
bit_rx = demapping(symb_rx_corected_down, Nbps, ModType); 

%% ====================== Generate Plots  =======================
bits_to_plot = min(params.timing.NumBits, 100 * Nbps); 
plotConstellation_Tx_Rx(ModOrder, ModType, symb_rx_down, symb_rx_corected_down);
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);
plotPSD_Tx_Rx(signal_tx_filtered, signal_rx_matched_filtered, Fs);
plotBasebandFrequencyResponse(signal_tx_filtered, signal_rx_matched_filtered, Fs);
plotFilterCharacteristics(g_rrc, Beta, Fs, OSF);

%% ====================== Robustness Analysis  =======================
EbN0_domain = 0:2:16;
num_iterations = 1000;

K_values = [1, 8, 16];
N_fixed = floor(length(symb_tx) * 0.1);
N_fixed = 20;
analyzeAvgWindow(params, pilot_position, EbN0_domain, num_iterations, ppm, K_values, N_fixed);

N_values = [10, 20, 40];
K_fixed = 8;
analyzePilotSize(params, pilot_position, EbN0_domain, num_iterations, ppm, N_values, K_fixed);

ppm_values = [0 2 20];
N = 40;
K = 8;
timing_offset_percent = 0.03;
cfoRobustness(params, pilot_position, EbN0_domain, num_iterations, ppm_values, N, K, timing_offset_percent);