function signal_out = applyFilter(signal_in, h_filter, NumFilterTaps)
    %   APPLYFILTER Applies an FIR filter using convolution and compensates for delay.
    %   Uses 'full' convolution and extracts the center part matching input length.
    %   Assumes h_filter is (approximately) linear phase.
    %
    %   Inputs:
    %       signal_in      - Input signal vector.
    %       h_filter       - Filter coefficients (impulse response).
    %       NumFilterTaps  - Length of h_filter.
    %
    %   Output:
    %       signal_out     - Filtered output signal (same length as input).
    

        % Perform 'full' convolution
        signal_filtered_full = conv(signal_in, h_filter, 'full');

        % Calculate group delay (approximate)
        delay = floor(NumFilterTaps / 2);

        % Extract the valid portion to compensate for delay and match input length
        signal_out = signal_filtered_full(delay + (1:length(signal_in)));
end