function out = downSampler(in, OSF)
    %   Downsamples the input signal by keeping every OSF-th sample, starting
    %   from the OSF-th sample (decimation).
    %
    %   Inputs:
    %       in  - Input row vector signal.
    %       OSF - Oversampling Factor (integer >= 1). Decimation factor.
    %
    %   Output:
    %       out - Downsampled output row vector (length = floor(length(in) / OSF)).

    % Select samples at indices OSF, 2*OSF, 3*OSF, ... up to the input length
    out = in(OSF:OSF:end);
end