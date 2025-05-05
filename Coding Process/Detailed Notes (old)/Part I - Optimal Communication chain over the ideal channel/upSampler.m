function out = upSampler(in, OSF)
    %UPSAMPLER Upsamples the input signal by inserting zeros.
    %   OUT = UPSAMPLER(IN, OSF) inserts OSF-1 zeros between each sample
    %   of the input row vector IN. This increases the sampling rate by a
    %   factor of OSF (Oversampling Factor).
    %
    %   Inputs:
    %       in  - Input row vector signal (e.g., sequence of symbols).
    %       OSF - Oversampling Factor (integer >= 1). The factor by which to
    %             increase the sampling rate. OSF=1 results in no change.
    %             Typically OSF >= 2 for practical pulse shaping.
    %
    %   Output:
    %       out - Output row vector signal, which is the upsampled version of IN.
    %             The length of OUT is length(IN) * OSF.
    %
    %   Example:
    %       in = [1, 2, 3];
    %       OSF = 3;
    %       out = upSampler(in, OSF); % out will be [0, 0, 1, 0, 0, 2, 0, 0, 3]
    %       (Note: This implementation places the sample at index OSF*i, so
    %        out would actually be [0 0 1 0 0 2 0 0 3] if N=3. For i=1, out(3)=in(1).
    %        For i=2, out(6)=in(2). For i=3, out(9)=in(3)).
    
        % Get the number of samples in the input signal.
        N = length(in);
    
        % Pre-allocate the output vector with zeros.
        % The output length will be N * OSF.
        out = zeros(1, OSF * N);
    
        % Iterate through each sample of the input signal.
        for i = 1:N
            % Place the i-th input sample at the corresponding position in the output vector.
            % The position is i * OSF, effectively inserting OSF-1 zeros before it
            % relative to the position (i-1)*OSF of the previous sample.
            % Example: If OSF=4, in(1) goes to out(4), in(2) goes to out(8), etc.
            out(OSF * i) = in(i);
        end
    end