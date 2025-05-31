%% ======= Step 3_v1 - Raw Implementation of Gardner Algorithm and basic Plots =========
%
% Implementationof Gardner Algorithm to overcome sampling time errors
%
%% ============================================================================================
clear; close all; clc;
addpath('../Part I - Optimal Communication chain over the ideal channel/p1_functions', 'p2_functions');



%% ========================================== Load Simulation Parameters  ==========================================
Nbps    = 2;
params  = initParameters(Nbps);
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

% ---- CFO and Sample Time Offset Parameters ----
Fc = 600e6;                                % Carrier frequency in Hz
delta_cfo_ppm   = 0.0 * 1e-6 * Fc;         % Frequency offset in Hz (0.08 ppm)
delta_omega     = 2 * pi * delta_cfo_ppm;  % Frequency offset in rad/s
phi_0           = 0;                       % Phase offset in rad
timing_offset_norm = 0.1;                  % Normalized timing offset (e.g., 0.9 of a symbol period)
initial_offset_samples = (timing_offset_norm * OSF); % Initial offset in samples introduced by circshift

% --- Gardner Algorithm Initialization ---
kappa = 0.08;                             % Gain for the Gardner algorithm.
current_base_sample_idx = ceil(OSF/2);
epsilon_hat = 0;                          % Initial estimate of the timing error (in samples).
symbols_corrected = zeros(Nsymb, 1);      % Array to store fixed output symbols
epsilon_k_history = zeros(Nsymb, 1);      % Array to store history of epsilon_hat



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
EbN0dB = 1000; % High SNR to observe Gardner clearly,
signal_tx_noisy = addAWGN(signal_tx, Eb, EbN0dB, OSF, SymRate);

% --- Introduce CFO, phase offset and sample time offset ---
num_samples_tx  = length(signal_tx_noisy);
time_vector     = (0 : num_samples_tx - 1).' * Ts;
offset_signal   = exp(1j * (delta_omega * time_vector + phi_0));
signal_tx_offset = signal_tx_noisy .* offset_signal;
signal_tx_offset = circshift(signal_tx_offset, initial_offset_samples);


% --- Receiver Chain ---
signal_rx   = applyFilter(signal_tx_offset, h_rrc, NumTaps);

% --- Gardner Loop ----
k=1;                                                                                    % Handle the first symbol separately
actual_sample_point_k1 = current_base_sample_idx + epsilon_hat;                         % Calculate initial FLAOTING POINT sampling instant for symbol 1

symbols_corrected(k) = interpolate_signal(signal_rx, actual_sample_point_k1);
epsilon_k_history(k) = epsilon_hat;                                                     % Store the initial timing error for symbol 1
y_prev_corr = symbols_corrected(k);                                                     % Set y_prev_corr for the next iteration (k=2), this is y[n-1]
current_base_sample_idx = current_base_sample_idx + OSF;                                % Advance the base sample index by OverSamplingFactor for the next symbol

for k = 2:Nsymb                                                                         % Loop for remaining symbols (n=2, 3, ... Nsymb)
    sample_point_k_minus_half = current_base_sample_idx - OSF/2 + epsilon_hat;
    sample_point_k = current_base_sample_idx + epsilon_hat;

    y_k_minus_half = interpolate_signal(signal_rx, sample_point_k_minus_half);          % Get mid-point sample y[n-0.5]
    y_k = interpolate_signal(signal_rx, sample_point_k);                     % Get current symbol sample y[n]

    err = real( y_k_minus_half * (conj(y_k) - conj(y_prev_corr)) );
    epsilon_hat = epsilon_hat - kappa * err;

    symbols_corrected(k) = y_k;                                     % Store the current symbol (sampled with epsilon_hat before this update)
    epsilon_k_history(k) = epsilon_hat;                             % Store the new timing error estimate

    y_prev_corr = y_k;                                              % Update y_prev_corr for the next iteration
    current_base_sample_idx = current_base_sample_idx + OSF;        % Advance base sample index for the next symbol
end


symb_rx = symbols_corrected;                                        % Assign corrected symbols to receiver output
bit_rx     = demapping(symb_rx, Nbps, ModType);                     % Demap symbols to bits
bit_rx     = bit_rx(:).';                                           % Ensure bit_rx is a row vector

%% ====================== Generate Plots  =======================
bits_to_plot = min(params.timing.NumBits, 100 * Nbps);
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, symb_rx);
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);
plotFilterCharacteristics(h_rrc, Beta, Fs, OSF);
plotPSD_Tx_Rx(signal_tx, signal_rx, Fs);

plotBasebandFrequencyResponse(signal_tx, signal_rx, Fs);

figure;
plot(1:Nsymb, epsilon_k_history);
title('Gardner Algorithm Timing Error Estimate (\epsilon_k)');
xlabel('Symbol Index (k)');
ylabel('Estimated Timing Error \epsilon_k (samples)');
grid on;


%% ====================== Interpolation Function =======================
function y = interpolate_signal(signal, float_index)
    N = length(signal);
    k = floor(float_index);
    mu = float_index - k;
    if k >= 1 && k < N
        y = (1 - mu) * signal(k) + mu * signal(k + 1);
    else
        y = 0; % Handle edge cases by zeroing out-of-bounds samples
    end
end