function params = initializeParameters()
    %INITIALIZEPARAMETERS Creates struct with all simulation parameters.
    %   Sets up modulation, timing, filter, sampling, and simulation settings.
    %   Calculates derived parameters and prints summary.
    %
    %   Outputs:
    %       params - Struct containing all parameters.

        % =====================================================================
        % == Modulation Parameters ==
        % =====================================================================
        params.modulation.Nbps = 4;                             % Bits per symbol (k)
        params.modulation.ModulationOrder = 2^params.modulation.Nbps; % M = 2^k
        % Determine Modulation Type based on Nbps (PAM for k=1 or odd, QAM for even k>=2)
        if params.modulation.Nbps == 1 || mod(params.modulation.Nbps, 2) ~= 0
            params.modulation.ModulationType = 'pam';
        else
            params.modulation.ModulationType = 'qam';
        end

        % =====================================================================
        % == Timing and Rate Parameters ==
        % =====================================================================
        params.timing.NumBits = params.modulation.Nbps * 100;   % Total data bits (multiple of Nbps)
        params.timing.SymbolRate = 5e6;                         % Symbol rate (Rs) [Hz]
        params.timing.SymbolPeriod = 1 / params.timing.SymbolRate; % Ts [s]
        params.timing.BitRate = params.timing.SymbolRate * params.modulation.Nbps; % Rb [bps]
        params.timing.BitPeriod = 1 / params.timing.BitRate;     % Tb [s]
        params.timing.StreamDuration = params.timing.NumBits / params.timing.BitRate; % Total time [s]

        % =====================================================================
        % == Filter Parameters ==
        % =====================================================================
        params.filter.RolloffFactor = 0.2;                      % RRC Roll-off factor (Beta)
        params.filter.NumFilterTaps = 101;                      % RRC Filter length (odd recommended)
        % Two-sided signal bandwidth: BW = Rs * (1 + Beta) [Hz]
        params.filter.SignalBandwidth = (1 + params.filter.RolloffFactor) * params.timing.SymbolRate;

        % =====================================================================
        % == Sampling Parameters ==
        % =====================================================================
        params.sampling.OversamplingFactor = 4;                 % Samples per symbol (OSF >= 2)
        params.sampling.SamplingFrequency = params.sampling.OversamplingFactor * params.timing.SymbolRate; % Fs [Hz]
        params.sampling.SamplePeriod = 1 / params.sampling.SamplingFrequency; % Tsamp [s]

        % =====================================================================
        % == BER Curve Simulation Parameters ==
        % =====================================================================
        params.simulation.EbN0_min_dB = 0;                      % Start Eb/N0 [dB]
        params.simulation.EbN0_max_dB = 15;                     % End Eb/N0 [dB]
        params.simulation.EbN0_step_dB = 1;                     % Step size Eb/N0 [dB]
        params.simulation.iterations_per_EbN0 = 50;             % Averaging iterations per point

        % =====================================================================
        % == Print Parameters to Command Window ==
        % =====================================================================
        fprintf('------------------------------------\n');
        fprintf('Simulation Parameters Initialized:\n');
        fprintf('------------------------------------\n');
        fprintf('  Modulation:\n');
        fprintf('    Bits per Symbol (Nbps):    %d\n', params.modulation.Nbps);
        fprintf('    Modulation Order (M):      %d\n', params.modulation.ModulationOrder);
        fprintf('    Modulation Type:           %s\n', upper(params.modulation.ModulationType));
        fprintf('  Timing & Rate:\n');
        fprintf('    Symbol Rate (Rs):          %.2e Symbols/sec\n', params.timing.SymbolRate);
        fprintf('    Symbol Period (Ts):        %.2e sec\n', params.timing.SymbolPeriod);
        fprintf('    Total Data Bits:           %d\n', params.timing.NumBits);
        fprintf('    Bit Rate (Rb):             %.2e bps\n', params.timing.BitRate);
        fprintf('    Total Signal Duration:     %.2e sec\n', params.timing.StreamDuration);
        fprintf('  Sampling:\n');
        fprintf('    Oversampling Factor (OSF): %d\n', params.sampling.OversamplingFactor);
        fprintf('    Sampling Frequency (Fs):   %.2e Hz\n', params.sampling.SamplingFrequency);
        fprintf('    Sample Period (Tsamp):     %.2e sec\n', params.sampling.SamplePeriod);
        fprintf('  Filtering (RRC):\n');
        fprintf('    Roll-off Factor (Beta):    %.2f\n', params.filter.RolloffFactor);
        fprintf('    Number of Taps:            %d\n', params.filter.NumFilterTaps);
        fprintf('    Signal Bandwidth (approx): %.2e Hz\n', params.filter.SignalBandwidth);
        fprintf('  BER Simulation:\n');
        fprintf('    Eb/N0 Range:               %.1f dB to %.1f dB\n', ...
                params.simulation.EbN0_min_dB, params.simulation.EbN0_max_dB);
        fprintf('    Eb/N0 Step:                %.1f dB\n', params.simulation.EbN0_step_dB);
        fprintf('    Iterations per Eb/N0:      %d\n', params.simulation.iterations_per_EbN0);
end