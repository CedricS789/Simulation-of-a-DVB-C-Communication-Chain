%% ======= Step 3_v1 - Raw Implementation of Gardner Algorithm =========

% Implementationof Gardner Algorithm to overcome sampling time errors
%
%% ============================================================================================

clear; close all; clc;
addpath('../Part I - Optimal Communication chain over the ideal channel/functions');
addpath('functions')

%% ========================================== Load Simulation Parameters  ==========================================
Nbps    = 2;
params  = initParameters_v2(Nbps);
NumBits = params.timing.NumBits;
ModType = params.modulation.ModulationType;
ModOrder= params.modulation.ModulationOrder;
OSF     = params.sampling.OversamplingFactor;
SymRate = params.timing.SymbolRate;
Tsymb   = params.timing.SymbolPeriod;
BitRate = params.timing.BitRate;
Fs      = params.sampling.SamplingFrequency;
Ts      = params.sampling.SamplePeriod;
Beta    = params.filter.RolloffFactor;
NumTaps = params.filter.NumFilterTaps;
iterations = params.simulation.iterations_per_EbN0;
displayParameters(params);

% ---- CFO Parameters ----
Fc = 600e6;                                 % Carrier frequency in Hz
delta_cfo_hz    = 0 * 1e-6 * Fc;            % Frequency offset in Hz (1 ppm)
delta_omega     = 2 * pi * delta_cfo_hz;    % Frequency offset in rad/s
phi_0           = 0;                        % Phase offset in rad

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

% --- Introduce CFO and phase offset ---
num_samples_tx  = length(signal_tx);                                % Number of samples in the transmitted signal
time_vector     = (0 : num_samples_tx - 1).' * Ts;                  % The TA insisted on this
offset_signal   = exp(1j * (delta_omega * time_vector + phi_0));    % Create the offset signal
signal_tx_offset   = signal_tx .* offset_signal;                    % Apply CFO to the transmitted signal
