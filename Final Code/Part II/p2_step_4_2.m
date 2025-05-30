%% =================== Step 4 - Frame and Frequency Acquisition ===================
%
% =================================================================================

clear; close all; clc;
addpath('../Part I - Optimal Communication chain over the ideal channel/p1_functions');
addpath('p2_functions');


%% ========================================== Load Simulation Parameters  ==========================================
Nbps    = 2;
params  = initParameters(Nbps);
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

% ---- CFO and sample time offset Parameters ----
Fc = 600e6;                                 % Carrier frequency in Hz
ppm = 2;
delta_cfo = ppm * 1e-6 * Fc;            % Frequency offset in Hz (1 ppm)
phi_0  = 0;                        % Phase offset in rad
sample_time_offset = 0;                   % Sample time offset

% ---- Data Segmentation Parameters ----
pilot_position = 50;
pilot_size = 20;
averaging_window = 8; % Number of symbols to average over for CFO estimation


%% ========================================== Communication Chain ==========================================
% --- Transmitter  ---
bit_tx      = randi([0, 1], 1, NumBits).';
symb_tx     = mapping(bit_tx, Nbps, ModType);
unuseful    = symb_tx(1 : pilot_position - 1);
pilot       = symb_tx(pilot_position : pilot_position + pilot_size - 1);
signal_tx   = upSampler(symb_tx, OSF).';
h_rrc       = rrcFilter(Beta, SymRate, OSF, NumTaps);
signal_tx_filtered   = applyFilter(signal_tx, h_rrc, NumTaps);
signalPower_tx  = mean(abs(signal_tx_filtered).^2);
Eb              = signalPower_tx / BitRate;

% -- Introduce Noise --
EbN0dB = 1000;
signal_tx_noisy = addAWGN(signal_tx_filtered, Eb, EbN0dB, OSF, SymRate);

% --- Introduce CFO, phase offset and sample time offset ---
time_vector = (0 : length(signal_tx) - 1).' * Ts;                  % The TA insisted on this
signal_tx_offset  = signal_tx_noisy .* exp(1j * (2 * pi * delta_cfo * time_vector + phi_0));   % Apply CFO to the transmitted signal


% --- Receiver Chain ---
signal_rx_matched_filtered  = applyFilter(signal_tx_offset, h_rrc, NumTaps);                % Apply Matched Filter
symb_rx    = downSampler(signal_rx_matched_filtered, OSF/2).';
[symb_rx, time_error] = gardner(symb_rx.');
[toa_hat, cfo_hat] = frameFreqAcquisition(pilot, symb_rx, averaging_window, Tsymb);
time_vector_symb = (0 : length(symb_rx) - 1).' * Tsymb;
symb_rx = symb_rx.' .* exp(-1j * (2 * pi * cfo_hat * time_vector_symb));  % Compensate for CFO and phase offset
bit_rx     = demapping(symb_rx, Nbps, ModType); 


%% ====================== Generate Plots  =======================
bits_to_plot = min(params.timing.NumBits, 100 * Nbps); 
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, symb_rx);
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);
plotFilterCharacteristics(h_rrc, Beta, Fs, OSF);
plotPSD_Tx_Rx(signal_tx_filtered, signal_rx_matched_filtered, Fs);
plotBasebandFrequencyResponse(signal_tx_filtered, signal_rx_matched_filtered, Fs);



%% ===================== Robustness Analysis ====================
averaging_window_vector = (1:1:20).';
toa_hat_vector = zeros(length(averaging_window_vector), 1);
cfo_hat_vector = zeros(length(averaging_window_vector), 1);
for i = 1:length(averaging_window_vector)
    [toa_hat_vector(i), cfo_hat_vector(i)] = frameFreqAcquisition(pilot, symb_rx, averaging_window_vector(i), Tsymb);
end

figure;
plot(averaging_window_vector, abs(toa_hat_vector - pilot_position));
xlabel('Averaging Window Size');
ylabel('abs(toa/_hat - pilot/_position)');
title('TOA Estimation Error vs. Averaging Window');
grid on;

figure;
plot(averaging_window_vector, abs(cfo_hat_vector - delta_cfo));
xlabel('Averaging Window Size');
ylabel('abs(.)');
title('TOA Estimation Error vs. Averaging Window');
grid on;


pilot_size_vector = (1:1:60).';
toa_hat_vector = zeros(length(pilot_size_vector), 1);
cfo_hat_vector = zeros(length(averaging_window_vector), 1);

K = 8;

for j = 1:length(pilot_size_vector)
    pilot = symb_tx(pilot_position : pilot_position + pilot_size_vector(j) - 1);
    [toa_hat_vector(j), cfo_hat_vector(j)] = frameFreqAcquisition(pilot, symb_rx, K, Tsymb);
end
figure;
plot(pilot_size_vector, abs(toa_hat_vector - pilot_position));
xlabel('Pilot Size (Number of Symbols)');
ylabel('abs(toa/_hat - pilot/_position)');
title('TOA Estimation Error vs. Pilot Size');
grid on;

figure;
plot(pilot_size_vector, abs(cfo_hat_vector - delta_cfo));
xlabel('Pilot Size (Number of Symbols)');
ylabel('abs(toa/_hat - pilot/_position)');
title('TOA Estimation Error vs. Pilot Size');
grid on;