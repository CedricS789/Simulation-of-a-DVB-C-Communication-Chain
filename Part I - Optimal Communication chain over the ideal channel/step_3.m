%%   =================== Step 3 - Simulation for BER Curve Generation ===================
%
%       Orchestrates the simulation over a range of Eb/N0 values to generate
%       a BER curve by simulating Tx -> AWGN Channel -> Rx repeatedly.
%
%   =======================================================================================

clear; close all; clc;

%% ========================================== Simulation Parameters  ==========================================
params = initializeParameters(); % Initialize fixed parameters

% --- Extract parameters needed ---
Nbps        = params.modulation.Nbps;
NumBits     = params.timing.NumBits;
ModType     = params.modulation.ModulationType;
ModOrder    = params.modulation.ModulationOrder;
OSF         = params.sampling.OversamplingFactor;
SymRate     = params.timing.SymbolRate;
Fs          = params.sampling.SamplingFrequency;
Beta        = params.filter.RolloffFactor;
NumTaps     = params.filter.NumFilterTaps;
BitRate     = params.timing.BitRate;

% --- BER Curve Parameters ---
EbN0_min_dB         = params.simulation.EbN0_min_dB;
EbN0_max_dB         = params.simulation.EbN0_max_dB;
EbN0_step_dB        = params.simulation.EbN0_step_dB;
iterations_per_EbN0 = params.simulation.iterations_per_EbN0; % Iterations for averaging BER

EbN0_domain_dB  = EbN0_min_dB:EbN0_step_dB:EbN0_max_dB; % Range of Eb/N0 values (dB)
num_EbN0_points = length(EbN0_domain_dB);               % Number of points on the curve

% --- Pre-allocate results array ---
ber_values = zeros(1, num_EbN0_points); % Stores simulated BER for each Eb/N0

fprintf('\n\n------------------------------------');
fprintf('\nBER Curve Simulation Setup:');
fprintf('\n------------------------------------');
fprintf('\n  Eb/N0 Range: %.1f dB to %.1f dB (Step: %.1f dB)', EbN0_min_dB, EbN0_max_dB, EbN0_step_dB);
fprintf('\n  Iterations per Eb/N0 point: %d', iterations_per_EbN0);
fprintf('\n  Total Bits per Eb/N0 point: %d', NumBits * iterations_per_EbN0);


%% ========================================== Transmitter Chain (Executed Once) ==========================================
fprintf('\n\n------------------------------------');
fprintf('\nTransmitter Processing (Once)');
fprintf('\n------------------------------------');

% -------- 1. Generate Source Data Bits --------
bit_tx = randi([0, 1], 1, NumBits);
fprintf('\n  Generated %d random bits.', NumBits);

% -------- 2. Symbol Mapping --------
symb_tx = mapping(bit_tx, Nbps, ModType);
fprintf('\n  Mapped bits to %d %d-%s symbols.', length(symb_tx), ModOrder, upper(ModType));

% -------- 3. Upsampling --------
symb_tx_up = upSampler(symb_tx, OSF).'; % Transpose for filter
fprintf('\n  Upsampled symbols by a factor of %d.', OSF);

% -------- 4. Pulse Shaping (Transmit Filtering) --------
h_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps); % Generates filter + optional plots
fprintf('\n  Generated RRC filter (beta=%.2f, %d taps).', Beta, NumTaps);
signal_tx = applyFilter(symb_tx_up, h_rrc, NumTaps);
fprintf('\n  Applied RRC pulse shaping filter.');

% -------- 5. Calculate Signal Metrics --------
signalPower_tx = mean(abs(signal_tx).^2);
Eb_tx = signalPower_tx / BitRate; % Energy per bit (Ptx / Rb) - constant for simulation
fprintf('\n  Calculated Tx Signal Power: %.2e', signalPower_tx);
fprintf('\n  Calculated Energy per Bit (Eb): %.2e', Eb_tx);
fprintf('\n  Transmitter Processing Complete.');


%% ========================================== Simulation Loop over Eb/N0 ==========================================
fprintf('\n\n------------------------------------');
fprintf('\nStarting Eb/N0 Loop...');
fprintf('\n------------------------------------');


% Outer Loop: Iterate through each specified Eb/N0 value
for idx_EbN0 = 1:num_EbN0_points
    EbN0dB = EbN0_domain_dB(idx_EbN0); % Current Eb/N0 in dB

    total_bit_errors = 0; % Reset error counter for this Eb/N0

    % --- Inner Loop for Averaging ---
    % Run multiple iterations for the same EbN0 to average out noise randomness
    for iter = 1:iterations_per_EbN0
        % -------- 1. AWGN Channel --------
        % Add noise based on current EbN0dB
        signal_tx_noisy = addAWGN(signal_tx, Eb_tx, EbN0dB, Fs);

        % -------- 2. Receiver Chain (Inlined) --------
        signal_rx_filtered = applyFilter(signal_tx_noisy, h_rrc, NumTaps); % Matched Filter
        symb_rx = downSampler(signal_rx_filtered, OSF).'; % Downsample
        bit_rx = demapping(symb_rx, Nbps, ModType); % Demap
        bit_rx = bit_rx(:).'; % Reshape to a vector

        % -------- 3. Calculate Errors for this iteration --------
        num_errors_iter = sum(bit_tx ~= bit_rx);
        total_bit_errors = total_bit_errors + num_errors_iter; % Accumulate errors
    end % End inner averaging loop

    % -------- Calculate Average BER for this EbN0 --------
    total_bits_simulated = NumBits * iterations_per_EbN0;

    if total_bits_simulated > 0
        ber_values(idx_EbN0) = total_bit_errors / total_bits_simulated; % Store average BER
    else
        ber_values(idx_EbN0) = NaN; % Avoid division by zero
    end

    fprintf('\n  Eb/N0 = %5.1f dB: Accumulated Errors = %6d, Average BER = %.3e', ...
            EbN0dB, total_bit_errors, ber_values(idx_EbN0));

end % End outer Eb/N0 loop

fprintf('\n\n------------------------------------');
fprintf('\nEb/N0 Loop Complete.');
fprintf('\n------------------------------------');


%% ========================================== Plotting Results ==========================================
fprintf('\n\n------------------------------------');
fprintf('\nPlotting BER Curve...');
fprintf('\n------------------------------------');

% Call the plotting function with results and parameters
plotBERCurve(EbN0_domain_dB, ber_values, ModType, ModOrder, OSF, Beta);

fprintf('\n\n------------------------------------');
fprintf('\nPlotting complete.');
fprintf('\n------------------------------------');