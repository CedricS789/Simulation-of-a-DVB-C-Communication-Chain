function signal_out = applyFilter(signal_in, h_filter, NumFilterTaps)
    signal_filtered_full = conv(signal_in, h_filter, 'full');
    delay = floor(NumFilterTaps / 2);
    signal_out = signal_filtered_full(delay + (1:length(signal_in)));
end