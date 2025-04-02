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
        params.simulation.EbN0_max_dB = 100;                     % End Eb/N0 [dB]
        params.simulation.EbN0_step_dB = 1;                     % Step size Eb/N0 [dB]
        params.simulation.iterations_per_EbN0 = 50;             % Averaging iterations per point

        % =====================================================================
        % == Print Parameters to Command Window ==
        % =====================================================================
        fprintf('------------------------------------');
        fprintf('\nSimulation Parameters Initialized:');
        fprintf('\n------------------------------------');
        fprintf('\n  Modulation:');
        fprintf('\n    Bits per Symbol (Nbps):    %d', params.modulation.Nbps);
        fprintf('\n    Modulation Order (M):      %d', params.modulation.ModulationOrder);
        fprintf('\n    Modulation Type:           %s', upper(params.modulation.ModulationType));
        fprintf('\n  Timing & Rate:');
        fprintf('\n    Symbol Rate (Rs):          %.2e Symbols/sec', params.timing.SymbolRate);
        fprintf('\n    Symbol Period (Ts):        %.2e sec', params.timing.SymbolPeriod);
        fprintf('\n    Total Data Bits:           %d', params.timing.NumBits);
        fprintf('\n    Bit Rate (Rb):             %.2e bps', params.timing.BitRate);
        fprintf('\n    Total Signal Duration:     %.2e sec', params.timing.StreamDuration);
        fprintf('\n  Sampling:');
        fprintf('\n    Oversampling Factor (OSF): %d', params.sampling.OversamplingFactor);
        fprintf('\n    Sampling Frequency (Fs):   %.2e Hz', params.sampling.SamplingFrequency);
        fprintf('\n    Sample Period (Tsamp):     %.2e sec', params.sampling.SamplePeriod);
        fprintf('\n  Filtering (RRC):');
        fprintf('\n    Roll-off Factor (Beta):    %.2f', params.filter.RolloffFactor);
        fprintf('\n    Number of Taps:            %d', params.filter.NumFilterTaps);
        fprintf('\n    Signal Bandwidth:          %.2e Hz', params.filter.SignalBandwidth);
        fprintf('\n  BER Simulation:');
        fprintf('\n    Eb/N0 Range:               %.1f dB to %.1f dB', ...
                params.simulation.EbN0_min_dB, params.simulation.EbN0_max_dB);
        fprintf('\n    Eb/N0 Step:                %.1f dB', params.simulation.EbN0_step_dB);
        fprintf('\n    Iterations per Eb/N0:      %d', params.simulation.iterations_per_EbN0);
end