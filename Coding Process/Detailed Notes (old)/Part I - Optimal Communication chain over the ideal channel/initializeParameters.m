function params = initializeParameters()
    %INITIALIZEPARAMETERS Creates and initializes a struct containing all simulation parameters.
    %   PARAMS = INITIALIZEPARAMETERS() sets up the default parameters for the
    %   communication system simulation, organizes them into nested structs
    %   (modulation, timing, filter, sampling, simulation), calculates derived
    %   parameters, and prints the configuration to the command window for review.
    %
    %   The parameters define aspects like the modulation scheme, data rates,
    %   pulse shaping filter characteristics, sampling rates, and the range
    %   of signal-to-noise ratios for BER simulation.
    %
    %   Outputs:
    %       params - A struct containing all the simulation parameters, organized
    %                into the following sub-structs:
    %                - params.modulation: Modulation type, order, bits per symbol.
    %                - params.timing: Symbol rate, bit rate, number of bits, durations.
    %                - params.filter: RRC filter properties (rolloff, taps, bandwidth).
    %                - params.sampling: Oversampling factor, sampling frequency/period.
    %                - params.simulation: Eb/N0 range and steps for BER curve generation.
    %
    %   Example Usage:
    %       simParams = initializeParameters();
    %       % Access parameters like: simParams.modulation.ModulationOrder
    %       % or simParams.timing.SymbolRate

        % =====================================================================
        % == Modulation Parameters ==
        % =====================================================================
        % Defines the digital modulation scheme used.
        params.modulation.Nbps = 4;                             % Number of bits per symbol (k = log2(M)).
                                                                % Examples: 1=BPSK, 2=QPSK, 3=8-PSK/8-PAM, 4=16-QAM/16-PAM, 6=64-QAM, 8=256-QAM.
        % Calculate Modulation Order (M), the number of distinct symbols.
        params.modulation.ModulationOrder = 2^params.modulation.Nbps; % M = 2^k

        % Determine Modulation Type ('pam' or 'qam') based on Nbps.
        % Convention used here: PAM for 1 bit (BPSK) or odd bits per symbol, QAM for even bits per symbol >= 2.
        if params.modulation.Nbps == 1 || mod(params.modulation.Nbps, 2) ~= 0
            params.modulation.ModulationType = 'pam'; % Pulse Amplitude Modulation (1D)
        else
            params.modulation.ModulationType = 'qam'; % Quadrature Amplitude Modulation (2D)
        end

        % =====================================================================
        % == Timing and Rate Parameters ==
        % =====================================================================
        % Defines symbol rates, bit rates, and simulation duration aspects.
        params.timing.NumBits = params.modulation.Nbps * 100;   % Total number of data bits to generate and transmit in the simulation.
                                                                % Chosen to be a multiple of Nbps for convenience.
        params.timing.SymbolRate = 5e6;                         % Symbol rate (Rs) in symbols per second (Hz). Also known as baud rate.
        params.timing.SymbolPeriod = 1 / params.timing.SymbolRate; % Duration of one symbol (Ts) in seconds.

        % Calculate the overall Bit Rate (Rb).
        params.timing.BitRate = params.timing.SymbolRate * params.modulation.Nbps; % Bits per second (Hz). Rb = Rs * log2(M).
        params.timing.BitPeriod = 1 / params.timing.BitRate;     % Duration of one bit (Tb) in seconds. Tb = Ts / log2(M).

        % Calculate the total duration of the bit stream transmission.
        params.timing.StreamDuration = params.timing.NumBits / params.timing.BitRate; % Total time in seconds. T_total = NumBits * Tb.
        % Alternative calculation: StreamDuration = (NumBits / Nbps) * SymbolPeriod

        % =====================================================================
        % == Filter Parameters ==
        % =====================================================================
        % Defines the properties of the pulse shaping and matched filter (RRC).
        params.filter.RolloffFactor = 0.2;                      % Roll-off factor (Beta, α) for the RRC filter (0 <= Beta <= 1).
                                                                % Controls excess bandwidth. Lower beta means sharper transition, longer impulse response.
        params.filter.NumFilterTaps = 101;                      % Number of taps (coefficients) for the RRC FIR filter.
                                                                % Should be odd for linear phase (symmetric impulse response).
                                                                % Must be large enough to capture the filter's response adequately.
        % Calculate the theoretical bandwidth occupied by the RRC filtered signal.
        % BW = Rs/2 * (1 + Beta), where Rs/2 is the minimum Nyquist bandwidth.
        params.filter.SignalBandwidth = params.timing.SymbolRate * (1 + params.filter.RolloffFactor) / 2 * 2; % Two-sided BW = Rs * (1 + Beta) Hz
        % Simpler calculation for double-sided bandwidth based on symbol rate
        params.filter.SignalBandwidth = (1 + params.filter.RolloffFactor) * params.timing.SymbolRate; % Actual occupied RF bandwidth (Hz)


        % =====================================================================
        % == Sampling Parameters ==
        % =====================================================================
        % Defines the sampling rates used in the simulation.
        params.sampling.OversamplingFactor = 4;                 % Oversampling Factor (OSF). Number of samples per symbol period.
                                                                % Must be >= 2 to satisfy Nyquist for the baseband signal bandwidth. Higher OSF provides better waveform representation.
        % Calculate the simulation's Sampling Frequency (Fs).
        params.sampling.SamplingFrequency = params.sampling.OversamplingFactor * params.timing.SymbolRate; % Samples per second (Hz). Fs = OSF * Rs.
        % Calculate the Sample Period (Tsamp).
        params.sampling.SamplePeriod = 1 / params.sampling.SamplingFrequency; % Duration of one sample in seconds. Tsamp = 1 / Fs = Ts / OSF.

        % =====================================================================
        % == BER Curve Simulation Parameters ==
        % =====================================================================
        % Defines the parameters for running the Bit Error Rate simulation loop.
        params.simulation.EbN0_min_dB = 0;                      % Starting value for the Eb/N0 range in dB.
        params.simulation.EbN0_max_dB = 15;                     % Ending value for the Eb/N0 range in dB.
        params.simulation.EbN0_step_dB = 1;                     % Step size for iterating through Eb/N0 values in dB.
        % Number of independent simulation runs (iterations) for each Eb/N0 point.
        % Averaging BER over multiple iterations reduces statistical variance and
        % provides a more reliable estimate, especially at low BERs.
        params.simulation.iterations_per_EbN0 = 50;

        % =====================================================================
        % == Print Parameters to Command Window ==
        % =====================================================================
        % Display the configured parameters for user verification.
        fprintf('------------------------------------\n');
        fprintf('Simulation Parameters Initialized:\n');
        fprintf('------------------------------------\n');

        % --- Modulation Section ---
        fprintf('  Modulation:\n');
        fprintf('    Bits per Symbol (Nbps):    %d\n', params.modulation.Nbps);
        fprintf('    Modulation Order (M):      %d\n', params.modulation.ModulationOrder);
        fprintf('    Modulation Type:           %s\n', upper(params.modulation.ModulationType));

        % --- Timing & Rate Section ---
        fprintf('  Timing & Rate:\n');
        fprintf('    Symbol Rate (Rs):          %.2e Symbols/sec\n', params.timing.SymbolRate);
        fprintf('    Symbol Period (Ts):        %.2e sec\n', params.timing.SymbolPeriod);
        fprintf('    Total Data Bits:           %d\n', params.timing.NumBits);
        fprintf('    Bit Rate (Rb):             %.2e bps\n', params.timing.BitRate); % Added BitRate printout
        fprintf('    Total Signal Duration:     %.2e sec\n', params.timing.StreamDuration);

        % --- Sampling Section ---
        fprintf('  Sampling:\n');
        fprintf('    Oversampling Factor (OSF): %d\n', params.sampling.OversamplingFactor);
        fprintf('    Sampling Frequency (Fs):   %.2e Hz\n', params.sampling.SamplingFrequency);
        fprintf('    Sample Period (Tsamp):     %.2e sec\n', params.sampling.SamplePeriod);

        % --- Filtering Section ---
        fprintf('  Filtering (RRC):\n');
        fprintf('    Roll-off Factor (Beta):    %.2f\n', params.filter.RolloffFactor);
        fprintf('    Number of Taps:            %d\n', params.filter.NumFilterTaps);
        fprintf('    Signal Bandwidth (approx): %.2e Hz\n', params.filter.SignalBandwidth); % Clarified label

        % --- BER Simulation Section ---
        fprintf('  BER Simulation:\n');
        fprintf('    Eb/N0 Range:               %.1f dB to %.1f dB\n', ...
                params.simulation.EbN0_min_dB, params.simulation.EbN0_max_dB);
        fprintf('    Eb/N0 Step:                %.1f dB\n', params.simulation.EbN0_step_dB); % Separated step printout
        fprintf('    Iterations per Eb/N0:      %d\n', params.simulation.iterations_per_EbN0);
    end