function out = upSampler(in, OSF)
    %UPSAMPLER Upsamples the input signal by inserting OSF-1 zeros between samples.
    %   Places original samples at indices k*OSF.
    %
    %   Inputs:
    %       in  - Input row vector signal.
    %       OSF - Oversampling Factor (integer >= 1).
    %
    %   Output:
    %       out - Upsampled output row vector (length = length(in) * OSF).

        N = length(in);
        out = zeros(1, OSF * N); % Pre-allocate output
        out(OSF:OSF:end) = in; % Place input samples at positions OSF, 2*OSF, ...
end