function params = coreParameters(Nbps_input)
    %INITIALIZEPARAMETERS Creates struct with all simulation parameters.
    %   Sets up modulation, timing, filter, sampling, and simulation settings.
    %   Calculates derived parameters and prints summary.
    %
    %   Outputs:
    %       params - Struct containing all parameters.

        % =====================================================================
        % == Modulation Parameters ==
        % =====================================================================
        params.modulation.Nbps = Nbps_input;                            % Bits per symbol (k)
        params.modulation.ModulationOrder = 2^params.modulation.Nbps;   % M = 2^k

        % Determine Modulation Type based on Nbps (PAM for k=1 or odd, QAM for even k>=2)
        if params.modulation.Nbps == 1 || mod(params.modulation.Nbps, 2) ~= 0
            params.modulation.ModulationType = 'pam';
        else
            params.modulation.ModulationType = 'qam';
        end

        % =====================================================================
        % == Timing and Rate Parameters ==
        % =====================================================================
        params.timing.NumBits = params.modulation.Nbps * 1e6;                           % Total data bits (multiple of Nbps)
        params.timing.SymbolRate = 5e6;                                                 % Symbol rate (Rs) [Hz]
        params.timing.SymbolPeriod = 1 / params.timing.SymbolRate;                      % Ts [s]
        params.timing.BitRate = params.timing.SymbolRate * params.modulation.Nbps;      % Rb [bps]
        params.timing.BitPeriod = 1 / params.timing.BitRate;                            % Tb [s]
        params.timing.StreamDuration = params.timing.NumBits / params.timing.BitRate;   % Total time [s]

        % =====================================================================
        % == Filter Parameters ==
        % =====================================================================
        params.filter.RolloffFactor = 0.2;                                                              % RRC Roll-off factor (Beta)
        params.filter.NumFilterTaps = 701;                                                              % RRC Filter length (odd recommended)
        params.filter.SignalBandwidth = (1 + params.filter.RolloffFactor) * params.timing.SymbolRate;   % Two-sided signal bandwidth: BW = Rs * (1 + Beta) [Hz]


        % =====================================================================
        % == Sampling Parameters ==
        % =====================================================================
        params.sampling.OversamplingFactor = 8;                                                             % Samples per symbol (OSF >= 2)
        params.sampling.SamplingFrequency = params.sampling.OversamplingFactor * params.timing.SymbolRate;  % Fs [Hz]
        params.sampling.SamplePeriod = 1 / params.sampling.SamplingFrequency; % Tsamp [s]

        % =====================================================================
        % == BER Curve Simulation Parameters ==
        % =====================================================================
        params.simulation.EbN0_min_dB = -5;                      % Start Eb/N0 [dB]
        params.simulation.EbN0_max_dB = 15;                      % End Eb/N0 [dB]
        params.simulation.EbN0_step_dB = 2;                      % Step size Eb/N0 [dB]
        params.simulation.iterations_per_EbN0 = 5;               % Averaging iterations per point

end