%%   =================== Step 2 v2 - Simulation with Fixed Noise Level (Noisy Channel) ===================
%
%       Purpose: This script extends the communication chain simulation to include an AWGN 
%       channel with a fixed Eb/N0. It evaluates the bit 
%       error rate (BER) averaged over multiple iterations of the noise realization (because the noise is random).
%
%       Context: Builds on Step 2 v1 by incorporating noise, maintaining the same 
%       symbol mapping and RRC filtering setup. The fixed Eb/N0 allows for an initial 
%       assessment of noise impact on BER, using a standardized energy per bit (Eb) definition 
%       based on signal power post-filtering. This iteration introduces performance evaluation 
%       while keeping the simulation scope limited to a single noise level.
%
%       Outputs: Produces BER statistics for the fixed Eb/N0, alongside plots of constellations, 
%       bit streams, and PSD, illustrating noise’s impact on signal integrity.
%
%   ======================================================================================

clear; close all; clc;
addpath('functions'); 



%% ========================== Initialization of Simulation Parameters ============================
Nbps = 1;
params = initParameters(Nbps); % Initialize parameters

% --- Extract parameters ---
NumBits     = params.timing.NumBits;
ModType     = params.modulation.ModulationType;
ModOrder    = params.modulation.ModulationOrder;
OSF         = params.sampling.OversamplingFactor;
SymRate     = params.timing.SymbolRate;
Fs          = params.sampling.SamplingFrequency;
Beta        = params.filter.RolloffFactor;
NumTaps     = params.filter.NumFilterTaps;
BitRate     = params.timing.BitRate;
displayParameters(params);

% --- Define Fixed Eb/N0 and Averaging Iterations ---
iterations = params.simulation.iterations_per_EbN0;    % Number of runs to average BER to have a good estimate of the BER
EbN0dB     = 30 ;                                       % The single Eb/N0 value (dB) for this test - ADJUST AS NEEDED

% ------------------------- Accurate BER Test Setup -------------------------
fprintf('\n\n========================================');
fprintf('\n         Accurate BER Test Setup        ');
fprintf('\n========================================');
fprintf('\n Modulation       : %d-%s', ModOrder, upper(ModType)); % Print modulation (e.g., 16-QAM)
fprintf('\n Target Eb/N0     : %.1f dB', EbN0dB);                 % Print target Eb/N0
fprintf('\n Bits per iter    : %d', NumBits);                     % Print bits per iteration
fprintf('\n Iterations       : %d', iterations);                  % Print number of iterations
fprintf('\n Total bits       : %d', NumBits * iterations);        % Print total bits simulated
fprintf('\n Symbol Rate      : %.3f MHz', SymRate/1e6);           % Print symbol rate in MHz
fprintf('\n Sampling Freq    : %.3f MHz', Fs/1e6);                % Print sampling frequency in MHz
fprintf('\n========================================\n'); 



%% ======================== Transmitter Chain (Executed Once) =====================
fprintf('\n\n========================================');
fprintf('\n         Transmitter Chain             ');       % Title for transmitter chain
fprintf('\n========================================');
fprintf('\n Generating RRC filter...');                     % Generate RRC filter...
h_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);             % Generate RRC filter

fprintf('\n Executing Transmitter Chain ...');             
bit_tx          = randi([0, 1], 1, NumBits);                % Generate random source bits
symb_tx         = mapping(bit_tx, Nbps, ModType);           % Map bits to complex symbols
symb_tx_up      = upSampler(symb_tx, OSF).';                % Upsample (transpose for filter)
signal_tx       = applyFilter(symb_tx_up, h_rrc, NumTaps);  % Apply RRC pulse shaping filter
signalPower_tx  = mean(abs(signal_tx).^2);                  % Average Tx signal power
Eb              = signalPower_tx / BitRate;                 % Energy per bit: Eb = P_avg / R_bit
fprintf('\n Measured Tx Signal Power (post-filter): %.2e', signalPower_tx);
fprintf('\n Transmitter Processing Complete.');
fprintf('\n========================================\n'); 




%% ================= Simulation Loop for Averaging at Fixed Eb/N0 ==================
fprintf('\n\n========================================');
fprintf('\n         Simulation Loop              ');
fprintf('\n========================================');
fprintf('\n Running %d iterations at Eb/N0 = %.1f dB...', iterations, EbN0dB);
total_bit_errors   = 0;                         % Accumulator for bit errors across iterations
total_bits_simulated = NumBits * iterations;    % Total bits for simulation

for iter = 1:iterations
    % --- AWGN Channel ---
    signal_tx_noisy = addAWGN(signal_tx, Eb, EbN0dB, OSF, SymRate);     % Add AWGN noise

    % --- Receiver Chain ---
    signal_rx  = applyFilter(signal_tx_noisy, h_rrc, NumTaps); % Matched filter
    symb_rx             = downSampler(signal_rx, OSF).';       % Downsample
    bit_rx              = demapping_v2(symb_rx, Nbps, ModType);         % Demap symbols
    bit_rx              = bit_rx(:).';                                  % Reshape to a vector

    % --- Calculate and Accumulate Errors ---
    num_errors_iter = sum(bit_tx ~= bit_rx);                            % Count bit errors for this iteration
    total_bit_errors = total_bit_errors + num_errors_iter; % Accumulate errors

    % Print iteration summary
    fprintf('\n Iteration %2d: Errors = %d', iter, num_errors_iter);
end
fprintf('\n Simulation Loop Complete.');
fprintf('\n========================================\n'); 



%% ==================== Performance Evaluation (Averaged Result) =====================
fprintf('\n\n========================================');
fprintf('\n       Performance Evaluation           ');
fprintf('\n========================================');
ber_avg = total_bit_errors / total_bits_simulated;                          % Average BER calculation
fprintf('\n Eb/N0       : %5.1f dB', EbN0dB);
fprintf('\n Eb          : %.2e', Eb);
fprintf('\n Total Bits  : %8d', total_bits_simulated);
fprintf('\n Total Errors: %6d', total_bit_errors);
fprintf('\n Average BER : %.3e', ber_avg);
fprintf('\n========================================\n'); 



%% ====================== Generate Plots (Based on Last Iteration) =======================
fprintf('\n\n========================================');
fprintf('\n         Generate Plots               ');
fprintf('\n========================================\n');
bits_to_plot = min(params.timing.NumBits, 100 * Nbps);                      % If NumBits is too large, plot only 100*Nbps bits
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, symb_rx);               % Plot Constellation
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);                          % Plot Bitstreams
plotPSD_Tx_Rx(signal_tx, signal_rx, Fs);                                    % Plot PSD
plotFilterCharacteristics(h_rrc, Beta, Fs, OSF);                            % Plot Filter Characteristics
plotBasebandFrequencyResponse(signal_tx, signal_rx, Fs);                     % Plot Baseband Frequency Response 

fprintf('\n========================================'); 
fprintf('\n           Plotting Complete            '); 
fprintf('\n========================================\n'); 
