%% ======= Step 3_v1 - Raw Implementation of Gardner Algorithm and Basic Plots =========
%
% Implementation of Gardner Algorithm to overcome sampling time errors
%
%% ============================================================================================

clear; close all; clc;
addpath('../Part I - Optimal Communication chain over the ideal channel/p1_functions');
addpath('p2_functions')

%% ========================================== Load Simulation Parameters  ==========================================
Nbps    = 4;
params  = initParameters_2(Nbps);
NumBits = params.timing.NumBits;
ModType = params.modulation.ModulationType;
ModOrder= params.modulation.ModulationOrder;
OSF     = params.sampling.OversamplingFactor;
SymRate = params.timing.SymbolRate;
Tsymb   = params.timing.SymbolPeriod;
Nsymb   = params.timing.NumSymbols;
BitRate = params.timing.BitRate;
Fs      = params.sampling.SamplingFrequency;
Ts      = params.sampling.SamplePeriod;
Beta    = params.filter.RolloffFactor;
NumTaps = params.filter.NumFilterTaps;
iterations = params.simulation.iterations_per_EbN0;

displayParameters(params);

% ---- CFO and Sample Time Offset Parameters ----
Fc = 600e6;                                 % Carrier frequency in Hz
delta_cfo_ppm   = 0.0 * 1e-6 * Fc;         % Frequency offset in Hz (0.08 ppm)
delta_omega     = 2 * pi * delta_cfo_ppm;   % Frequency offset in rad/s
phi_0           = 0;                        % Phase offset in rad
timing_offset_norm = 0.9;                   % Normalized timing offset (10% of symbol period)
initial_offset_samples = round(timing_offset_norm * OSF); % Initial offset in samples (1 sample)

% --- Gardner Algorithm Initialization ---
kappa = 0.01;                               % Gain for the Gardner algorithm
epsilon_hat = 0;                            % Initial timing offset estimate
y_n_minus_1 = complex(0, 0);                % Previous symbol sample initialized to zero
output_symbols = zeros(Nsymb, 1);           % Pre-allocate corrected symbols
epsilon_hat_history = zeros(Nsymb, 1);      % Track timing estimate convergence

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

% --- Introduce CFO, Phase Offset, and Sample Time Offset ---
num_samples_tx  = length(signal_tx_noisy);
time_vector     = (0 : num_samples_tx - 1).' * Ts;
offset_signal   = exp(1j * (delta_omega * time_vector + phi_0));
signal_tx_offset = signal_tx_noisy .* offset_signal;
signal_tx_offset = circshift(signal_tx_offset, initial_offset_samples); % Apply timing offset (1 sample delay)

% --- Receiver Chain ---
signal_rx  = applyFilter(signal_tx_offset, h_rrc, NumTaps);

% --- Gardner Algorithm for Timing Recovery ---
for n = 1 : Nsymb
    % Fractional index for current symbol y[n]
    frac_index_n = (n + epsilon_hat) * OSF;
    y_n = interpolate_signal(signal_rx, frac_index_n);
    
    if n >= 2
        % Fractional index for midway sample y[n-0.5]
        frac_index_mid = (n - 0.5 + epsilon_hat) * OSF;
        y_mid = interpolate_signal(signal_rx, frac_index_mid);
        
        % Timing error computation
        e = real(y_mid * (conj(y_n) - conj(y_n_minus_1)));
    else
        e = 0; % No error for first symbol (no y[n-1])
    end
    
    % Update timing estimate
    epsilon_hat = epsilon_hat - kappa * e;
    epsilon_hat_history(n) = epsilon_hat;
    
    % Store corrected symbol
    output_symbols(n) = y_n;
    
    % Update previous symbol
    y_n_minus_1 = y_n;
end

symb_rx = output_symbols; % Use Gardner-corrected symbols
bit_rx  = demapping(symb_rx, Nbps, ModType);
bit_rx  = bit_rx(:).';

%% ====================== Generate Plots  =======================
bits_to_plot = min(params.timing.NumBits, 100 * Nbps);
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, symb_rx);
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);
plotFilterCharacteristics(h_rrc, Beta, Fs, OSF);
plotPSD_Tx_Rx(signal_tx, signal_rx, Fs);
plotBasebandFrequencyResponse(signal_tx, signal_rx, Fs);

% --- Additional Plot for Gardner Convergence ---
figure;
plot(epsilon_hat_history, 'b-', 'LineWidth', 1.5);
title('Gardner Timing Estimate Convergence');
xlabel('Symbol Index (n)');
ylabel('Estimated Timing Offset (\epsilon[n])');
grid on;

%% ====================== Interpolation Function =======================
function y = interpolate_signal(signal, frac_index)
    N = length(signal);
    k = floor(frac_index);
    mu = frac_index - k;
    if k >= 1 && k < N
        y = (1 - mu) * signal(k) + mu * signal(k + 1);
    else
        y = 0; % Handle edge cases by zeroing out-of-bounds samples
    end
end