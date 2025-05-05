function signal_out = applyFilter(signal_in, h_filter, NumFilterTaps)
    %APPLYFILTER Applies an FIR filter to a signal using convolution and handles delay.
    %   SIGNAL_OUT = APPLYFILTER(SIGNAL_IN, H_FILTER, NUMFILTERTAPS) filters
    %   the input signal SIGNAL_IN with the FIR filter defined by the impulse
    %   response H_FILTER. It uses 'full' convolution and then extracts the
    %   central part of the result to maintain the same signal length as the
    %   input and compensate for the filter's group delay.
    %
    %   This method assumes H_FILTER represents a causal, linear-phase FIR
    %   filter (or approximates one), where the group delay is roughly half
    %   the filter length.
    %
    %   Inputs:
    %       signal_in      - Input signal vector (row or column).
    %       h_filter       - Filter coefficients (impulse response) as a vector.
    %       NumFilterTaps  - The number of taps (length) of the filter h_filter.
    %                        This must match length(h_filter). Passed explicitly
    %                        to avoid recalculating length inside.
    %
    %   Output:
    %       signal_out     - Output signal vector, filtered version of signal_in,
    %                        with the same dimensions and length as signal_in.
    %
    %   Example:
    %       signal = randn(1, 100);
    %       filter_coeffs = ones(1, 11)/11; % Simple moving average filter
    %       num_taps = 11;
    %       filtered_signal = applyFilter(signal, filter_coeffs, num_taps);
    %       % filtered_signal will have length 100.
    
        % Perform convolution between the input signal and the filter's impulse response.
        % The 'full' option computes the complete convolution, resulting in an output
        % vector of length length(signal_in) + length(h_filter) - 1.
        % This includes transient effects at the beginning and end.
        signal_filtered_full = conv(signal_in, h_filter, 'full');
    
        % Calculate the approximate group delay of the filter.
        % For a linear-phase FIR filter with NumFilterTaps taps, the delay is
        % (NumFilterTaps - 1) / 2 samples.
        % Using floor(NumFilterTaps / 2) provides the correct integer index offset
        % needed to center the output for both odd and even NumFilterTaps,
        % although odd is typical for symmetric filters.
        delay = floor(NumFilterTaps / 2);
    
        % Extract the relevant portion of the full convolution result.
        % We want the output signal to align in time with the input signal,
        % effectively removing the filter's delay. We skip the initial 'delay'
        % samples (transient part) and take the next 'length(signal_in)' samples.
        % The indexing `delay + (1:length(signal_in))` creates indices from
        % `delay + 1` to `delay + length(signal_in)`.
        signal_out = signal_filtered_full(delay + (1:length(signal_in)));
    end