%% ====================================== p2_step_3 - Gardner  ======================================
close all; clear; clc;
addpath('../Part I/p1_functions');
addpath('p2_functions')


%% =================== Load Simulation Parameters  ===================
Nbps = 2;
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


%% =================== Communication Chain ===================
% ---- CFO and Sample Time Offset Parameters ----
Fc = 600e6;                                                    % Carrier frequency in Hz
timing_offset_percent = 0.05;                                   % Normalized timing offset (% of symbol period)
initial_offset_samples = round(timing_offset_percent * OSF);   % Initial offset in samples (1 sample)

% --- Transmitter  ---
bit_tx = randi([0, 1], 1, NumBits).';
symb_tx = mapping(bit_tx, Nbps, ModType);
symb_tx_up = upSampler(symb_tx, OSF).';
g_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);
signal_tx_filtered = applyFilter(symb_tx_up, g_rrc, NumTaps);
signalPower_tx = mean(abs(signal_tx_filtered).^2);
Eb = signalPower_tx / BitRate;

time_vector = (0 : length(signal_tx_filtered) - 1).' * Ts;
time_vector_symb = (0 : length(symb_tx) - 1).' * Tsymb;

% -- Introduce Noise --
EbN0dB = 1e10;
signal_tx_noisy = addAWGN(signal_tx_filtered, Eb, EbN0dB, OSF, SymRate);

% --- Introduce CFO, Phase Offset, and Sample Time Offset ---
signal_tx_distorted = circshift(signal_tx_noisy, initial_offset_samples);

% --- Receiver Chain ---
kappa = 0.001;
signal_rx_matched_filtered  = applyFilter(signal_tx_distorted, g_rrc, NumTaps);
symb_rx_down = downSampler(signal_rx_matched_filtered, OSF);
[symb_rx_corected_down, time_shift_errors] = gardner(signal_rx_matched_filtered, kappa, OSF);
bit_rx_corrected = demapping(symb_rx_corected_down, Nbps, ModType);

%% =================== Generate Plots  ===================
plotConstellation_Tx_Rx(ModOrder, ModType, symb_rx_down, symb_rx_corected_down);
bits_to_plot = min(params.timing.NumBits, 100 * Nbps); 
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx_up, symb_rx_down);
plotBitstream_Tx_Rx(bit_tx, bit_rx_corrected, bits_to_plot);
plotPSD_Tx_Rx(signal_tx_filtered, signal_rx_matched_filtered, Fs);
plotBasebandFrequencyResponse(signal_tx_filtered, signal_rx_matched_filtered, Fs);
plotFilterCharacteristics(g_rrc, Beta, Fs, OSF);