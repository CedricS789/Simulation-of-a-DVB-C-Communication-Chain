function displayParameters(params)
% displayParameters prints the simulation parameters in a structured format.
%
%   displayParameters(params)
%
% Inputs:
%   params - Structure containing simulation parameters (from initParameters)

    % Print header
    fprintf('========================================\n');
    fprintf('Simulation Parameters Overview:\n');
    fprintf('========================================');
    
    % Modulation Parameters
    fprintf('\nModulation Parameters:\n');
    fprintf('  Bits per symbol      : %d\n', params.modulation.Nbps);
    fprintf('  Modulation order     : %d\n', params.modulation.ModulationOrder);
    fprintf('  Modulation type      : %s\n', upper(params.modulation.ModulationType));
    
    % Timing and Rate Parameters
    fprintf('\nTiming and Rate Parameters:\n');
    fprintf('  Number of Bits sent  : %d\n', params.timing.NumBits);
    fprintf('  Symbol rate          : %.2f MHz\n', params.timing.SymbolRate/1e6);
    fprintf('  Symbol period        : %.2e s\n', params.timing.SymbolPeriod);
    fprintf('  Bit rate             : %.2f bps\n', params.timing.BitRate);
    fprintf('  Bit period           : %.2e s\n', params.timing.BitPeriod);
    fprintf('  Stream duration      : %.2f s\n', params.timing.StreamDuration);
    
    % Filter Parameters
    fprintf('\nFilter Parameters:\n');
    fprintf('  Rolloff factor       : %.2f\n', params.filter.RolloffFactor);
    fprintf('  Number of filter taps: %d\n', params.filter.NumFilterTaps);
    fprintf('  Signal bandwidth     : %.2f MHz\n', params.filter.SignalBandwidth/1e6);
    
    % Sampling Parameters
    fprintf('\nSampling Parameters:\n');
    fprintf('  Oversampling factor  : %d\n', params.sampling.OversamplingFactor);
    fprintf('  Sampling frequency   : %.2f MHz\n', params.sampling.SamplingFrequency/1e6);
    fprintf('  Sample period        : %.2e s\n', params.sampling.SamplePeriod);
    
    % BER Curve Simulation Parameters
    fprintf('\nBER Curve Simulation Parameters:\n');
    fprintf('  Eb/N0 range          : [%g dB to %g dB] with step = %.1f dB\n', ...
            params.simulation.EbN0_min_dB, params.simulation.EbN0_max_dB, params.simulation.EbN0_step_dB);
    fprintf('  Eb/N0 values         : ');
    fprintf('%g ', params.simulation.EbN0_domain_dB);
    fprintf(' dB\n');
    fprintf('  Iterations per Eb/N0 : %d', params.simulation.iterations_per_EbN0);
    fprintf('\n========================================\n'); 

end
