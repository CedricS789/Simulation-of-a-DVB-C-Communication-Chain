function out = downSampler(in, OSF)
    %   Downsamples the input signal by keeping every OSF-th sample.
    %   OUT = DOWNSAMPLER(IN, OSF) selects every OSF-th sample from the input
    %   row vector IN, starting from the OSF-th sample. This reduces the
    %   sampling rate by a factor of OSF (Oversampling Factor), effectively
    %   performing decimation without an anti-aliasing filter.
    %
    %   Inputs:
    %       in  - Input row vector signal (e.g., received signal after matched filter).
    %       OSF - Oversampling Factor (integer >= 1). The factor by which to
    %             decrease the sampling rate. OSF=1 results in no change.
    %
    %   Output:
    %       out - Output row vector signal, which is the downsampled version of IN.
    %             The length of OUT is floor(length(IN) / OSF).
    %
    %   Example:
    %       in = [1, 2, 3, 4, 5, 6, 7, 8, 9];
    %       OSF = 3;
    %       out = downSampler(in, OSF); % out will be [3, 6, 9]
    %       (Selects samples at indices 3, 6, 9).
    %
    %   Note: This specific implementation selects samples at indices OSF, 2*OSF, 3*OSF, ...
    %         This corresponds to a specific sampling phase. Other downsampling phases
    %         (e.g., selecting samples 1, 1+OSF, 1+2*OSF, ...) might be needed
    %         depending on synchronization requirements. Assumes input signal 'in'
    %         has appropriate length related to OSF.
    
        % Calculate the number of output samples.
        % Use floor to handle cases where the input length is not perfectly divisible by OSF.
        N = floor(length(in) / OSF);
    
        % Pre-allocate the output vector with zeros.
        out = zeros(1, N);
    
        % Iterate through the indices of the output samples.
        for i = 1:N
            % Select the sample from the input vector at the index i * OSF.
            % This effectively keeps one sample for every OSF input samples.
            out(i) = in(i * OSF);
        end
    end