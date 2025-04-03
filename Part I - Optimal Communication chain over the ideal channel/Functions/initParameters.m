function params = initParameters(Nbps_input)
    %   Creates struct with all simulation parameters based on the desired modulation scheme (Number of bits per symbol).
    %   Sets up modulation, timing, filter, sampling, and simulation settings.
    %   Calculates derived parameters such as symbol rate, bit rate, and bandwidth etc.
    %   
    %   Inputs:
    %       Nbps_input - Number of bits per symbol (k) for modulation.
    %       Nbps_input = 1 for PAM (k=1), 2 for QPSK (k=2), 3 for 8-PSK (k=3), etc.
    %   Outputs:
    %       params - Struct containing all parameters for the modulation scheme (Nbps_input).

        % =====================================================================
        % == Modulation Parameters ==
        % =====================================================================
        params.modulation.Nbps = Nbps_input;                            % Bits per symbol (k)
        params.modulation.ModulationOrder = 2^params.modulation.Nbps;   % M = 2^k

        % Determine Modulation Type based on Nbps (PAM for k=1 or odd, QAM for even k>=2)
        if params.modulation.Nbps == 1 || mod(params.modulation.Nbps, 2) ~= 0       % if Npbs is odd then PAM
            params.modulation.ModulationType = 'pam';
        else
            params.modulation.ModulationType = 'qam';
        end

        % =====================================================================
        % == Timing and Rate Parameters ==
        % =====================================================================
        params.timing.NumBits = params.modulation.Nbps * 1e5;                           % Total data bits (multiple of Nbps)
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
        params.sampling.SamplePeriod = 1 / params.sampling.SamplingFrequency;                               % Tsamp [s]

        % =====================================================================
        % == BER Curve Simulation Parameters ==
        % =====================================================================
        params.simulation.EbN0_min_dB = -5;                      % Start Eb/N0 [dB]
        params.simulation.EbN0_max_dB = 15;                      % End Eb/N0 [dB]
        params.simulation.EbN0_step_dB = 1;                      % Step size Eb/N0 [dB]
        params.simulation.iterations_per_EbN0 = 1;               % Averaging iterations per point
        
        EbN0_min_dB                         = params.simulation.EbN0_min_dB;                % Minimum Eb/N0 value in dB
        EbN0_max_dB                         = params.simulation.EbN0_max_dB;                % Maximum Eb/N0 value in dB
        EbN0_step_dB                        = params.simulation.EbN0_step_dB;               % Step size for Eb/N0 sweep in dB
        params.simulation.EbN0_domain_dB    = EbN0_min_dB:EbN0_step_dB:EbN0_max_dB;         % Range of Eb/N0 values to simulate (dB)

end