%% ======= Step 3_v1 - Raw Implementation of Gardner Algorithm and basic Plots =========
%
% Implementationof Gardner Algorithm to overcome sampling time errors
%
%% ============================================================================================


clear; close all; clc;
addpath('../Part I - Optimal Communication chain over the ideal channel/p1_functions');
addpath('p2_functions')

%% ========================================== Load Simulation Parameters  ==========================================
Nbps    = 2;
params  = initParameters_v2(Nbps);
NumBits = params.timing.NumBits;
ModType = params.modulation.ModulationType;
ModOrder= params.modulation.ModulationOrder;
OSF     = params.sampling.OversamplingFactor;
SymRate = params.timing.SymbolRate;
Tsymb   = params.timing.SymbolPeriod;
Nsymb = params.timing.NumSymbols;
BitRate = params.timing.BitRate;
Fs      = params.sampling.SamplingFrequency;
Ts      = params.sampling.SamplePeriod;
Beta    = params.filter.RolloffFactor;
NumTaps = params.filter.NumFilterTaps;
iterations = params.simulation.iterations_per_EbN0;

displayParameters(params);

% ---- CFO and sample time offset Parameters ----
Fc = 600e6;                                 % Carrier frequency in Hz
delta_cfo_ppm   = 0.08 * 1e-6 * Fc;         % Frequency offset in Hz (1 ppm)
delta_omega     = 2 * pi * delta_cfo_ppm;   % Frequency offset in rad/s
phi_0           = 0;                        % Phase offset in rad

% --- Gardner Algorithm Initialization ---
kappa = 0.01;                               % Gain for the Gardner algorithm. This is the crucial tuning parameter. It controls how quickly the loop adapts. Too large, and it might overshoot or become unstable; too small, and it converges slowly.
timing_offset_norm = 0.1;           % Normalized timing offset (0.1 = 10% of the symbol period)

%% ========================================== Communication Chain ==========================================
% --- Transmitter  ---
bit_tx      = randi([0, 1], 1, NumBits).';
symb_tx     = mapping(bit_tx, Nbps, ModType);
symb_tx_up  = upSampler(symb_tx, OSF).';
h_rrc       = rrcFilter(Beta, SymRate, OSF, NumTaps);
signal_tx   = applyFilter(symb_tx_up, h_rrc, NumTaps);
signalPower_tx  = mean(abs(signal_tx).^2);
Eb              = signalPower_tx / BitRate;

% -- Introduce Noise --
EbN0dB = 1000;
signal_tx_noisy = addAWGN(signal_tx, Eb, EbN0dB, OSF, SymRate);

% --- Introduce CFO, phase offset and sample time offset ---
num_samples_tx  = length(signal_tx_noisy);                          % Number of samples in the transmitted signal
time_vector     = (0 : num_samples_tx - 1).' * Ts;                  % The TA insisted on this
offset_signal   = exp(1j * (delta_omega * time_vector + phi_0));    % Create the offset signal
signal_tx_offset = signal_tx_noisy .* offset_signal;                % Apply CFO to the transmitted signal
signal_tx_offset = circshift(signal_tx_offset, timing_offset_norm); % Represents th e time shift at receiver (1%)


% --- Receiver Chain ---
signal_rx  = applyFilter(signal_tx_offset, h_rrc, NumTaps);
symb_rx    = downSampler(signal_rx, OSF).';
bit_rx     = demapping(symb_rx, Nbps, ModType); 
bit_rx     = bit_rx(:).';  


%% ====================== Generate Plots  =======================
bits_to_plot = min(params.timing.NumBits, 100 * Nbps); 
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, symb_rx);
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);
plotFilterCharacteristics(h_rrc, Beta, Fs, OSF);
plotPSD_Tx_Rx(signal_tx, signal_rx, Fs);
plotBasebandFrequencyResponse(signal_tx, signal_rx, Fs);

