function hFig = plotFilterCharacteristics(h_rrc, Beta, Fs, OSF)
    %PLOTFILTERCHARACTERISTICS Visualizes RRC and combined RC filter properties.
    %   Includes RRC impulse/frequency response and combined RC pulse ISI check.
    %
    %   Inputs:
    %       h_rrc             - RRC filter coefficients (impulse response).
    %       Beta              - Roll-off factor (for annotation).
    %       Fs                - Sampling frequency (Hz) (for axis scaling).
    %       OSF               - Oversampling factor (for ISI check).
    %
    %   Output:
    %       hFig              - Handle to the created figure object.

    % =====================================================================
    % == Plotting Parameters (Locally Defined) ==
    % =====================================================================
    filterColor = [0.4660, 0.6740, 0.1880]; % Greenish
    linkColor   = [0.4940, 0.1840, 0.5560]; % Purple
    rxColor     = [0.8500, 0.3250, 0.0980]; % Reddish
    lineWidth = 1.5;
    baseMarkerSize = 6;
    rrcImpulseStemStyle = 'filled';
    rrcImpulseMarkerSize= 4;
    isiStemMarkerShape  = 'o';
    isiStemLineStyle    = 'none';
    figureName          = 'Filter Analysis (RRC & RC)';
    figurePosition      = [100, 100, 800, 700];
    useFigureNumberTitle= 'off';
    NFFT                = 2^12;             % FFT points for freq response calc
    NUM_SYMBOLS_SIDE_ISI= 5;                % Symbol periods +/- center for ISI plot
    timeUnitsFactor     = 1e6;              % s to us
    timeUnitsLabel      = '\mus';
    freqUnitsFactor     = 1e6;              % Hz to MHz
    freqUnitsLabel      = 'MHz';
    subplot1TitleFormat = 'RRC Filter Impulse Response (β=%.2f, Taps=%d)';
    subplot1XLabel      = sprintf('Time (%s)', timeUnitsLabel);
    subplot1YLabel      = 'Amplitude';
    subplot2Title       = 'RRC Filter Frequency Response (Magnitude)';
    subplot2XLabel      = sprintf('Frequency (%s)', freqUnitsLabel);
    subplot2YLabel      = 'Magnitude (dB)';
    freqResponseMinDb   = -80;              % Y-axis floor (dB) for freq response
    freqResponseMaxDb   = 10;               % Y-axis ceiling (dB) for freq response
    subplot3Title       = 'Combined Raised Cosine Pulse Response & ISI Check';
    subplot3XLabel      = 'Time (Normalized to Symbol Periods, t/T_{symb})';
    subplot3YLabel      = 'Amplitude';
    isiPlotXlimPadding  = 0.5;              % Padding (symbol periods) for ISI plot x-axis
    isiLegendLocation   = 'best';
    isiTextContent      = ' Samples at k*T_{symb} ';
    isiTextYPosFactor   = 0.9;              % Relative Y position for ISI text
    isiTextHAlign       = 'left';
    isiTextVAlign       = 'top';
    isiTextBGColor      = 'w';
    isiTextEdgeColor    = 'k';
    isiTextFontSize     = 8;
    % =====================================================================

    % =====================================================================
    % == Derived Parameters & Constants ==
    % =====================================================================
    NumFilterTaps = length(h_rrc);
    Ts = 1/Fs; % Sample period
    Tsymb = Ts * OSF; % Symbol period

    % =====================================================================
    % == Plotting Code ==
    % =====================================================================
    hFig = figure('Name', figureName, ...
                  'NumberTitle', useFigureNumberTitle, ...
                  'Position', figurePosition);

    % --- Subplot 1: RRC Impulse Response (Time Domain) ---
    ax1 = subplot(3, 1, 1, 'Parent', hFig);
    time_axis_filter_sec = (-(NumFilterTaps-1)/2 : (NumFilterTaps-1)/2) * Ts;
    time_axis_filter_units = time_axis_filter_sec * timeUnitsFactor;
    stem(ax1, time_axis_filter_units, h_rrc, ...
        rrcImpulseStemStyle, 'Color', filterColor, 'MarkerSize', rrcImpulseMarkerSize);
    grid(ax1, 'on'); box(ax1, 'on');
    xlabel(ax1, subplot1XLabel);
    ylabel(ax1, subplot1YLabel);
    title(ax1, sprintf(subplot1TitleFormat, Beta, NumFilterTaps));
    xlim(ax1, time_axis_filter_units([1 end]));

    % --- Subplot 2: RRC Frequency Response (Magnitude) ---
    ax2 = subplot(3, 1, 2, 'Parent', hFig);
    [H_full, ~] = freqz(h_rrc, 1, NFFT, 'whole', Fs);
    H_shifted = fftshift(H_full);
    F_shifted_hz = (-NFFT/2 : NFFT/2 - 1) * (Fs / NFFT);
    F_shifted_units = F_shifted_hz / freqUnitsFactor;
    plot(ax2, F_shifted_units, 20*log10(abs(H_shifted) + eps), ... % Use eps for stability
        'Color', filterColor, 'LineWidth', lineWidth);
    grid(ax2, 'on'); box(ax2, 'on');
    xlabel(ax2, subplot2XLabel);
    ylabel(ax2, subplot2YLabel);
    title(ax2, subplot2Title);
    xlim(ax2, [-Fs/2, Fs/2] / freqUnitsFactor);
    ylim(ax2, [freqResponseMinDb, freqResponseMaxDb]);

    % --- Subplot 3: Combined RC Response & ISI Check ---
    ax3 = subplot(3, 1, 3, 'Parent', hFig);
    h_rc = conv(h_rrc, h_rrc); % Combined Tx+Rx filter response
    NumRCTaps = length(h_rc);
    rc_center_index = floor(NumRCTaps / 2) + 1;
    time_axis_rc_sec = (-(NumRCTaps-1)/2 : (NumRCTaps-1)/2) * Ts;
    time_axis_rc_norm = time_axis_rc_sec / Tsymb; % Normalize time axis
    plot(ax3, time_axis_rc_norm, h_rc, 'Color', linkColor, 'LineWidth', lineWidth);
    hold(ax3, 'on');

    % --- ISI Check: Mark samples at integer multiples of Tsymb ---
    symbol_indices_offset = -NUM_SYMBOLS_SIDE_ISI:1:NUM_SYMBOLS_SIDE_ISI;
    % Calculate sample indices corresponding to symbol instants (k * Tsymb)
    symbol_sample_indices = rc_center_index + symbol_indices_offset * OSF;
    valid_indices_mask = (symbol_sample_indices >= 1) & (symbol_sample_indices <= NumRCTaps);
    symbol_sample_indices = symbol_sample_indices(valid_indices_mask);
    symbol_times_norm = time_axis_rc_norm(symbol_sample_indices);
    symbol_values = h_rc(symbol_sample_indices);
    % Plot stems at symbol sampling instants
    stem(ax3, symbol_times_norm, symbol_values, ...
         'Color', rxColor, 'Marker', isiStemMarkerShape, ...
         'MarkerFaceColor', rxColor, 'MarkerSize', baseMarkerSize, ...
         'LineStyle', isiStemLineStyle);
    legend(ax3, 'RC Pulse Shape', sprintf('Samples at k*T_{symb} (OSF=%d)', OSF), ...
           'Location', isiLegendLocation);
    grid(ax3, 'on'); box(ax3, 'on');
    xlabel(ax3, subplot3XLabel);
    ylabel(ax3, subplot3YLabel);
    title(ax3, subplot3Title);
    xlim(ax3, [-NUM_SYMBOLS_SIDE_ISI - isiPlotXlimPadding, NUM_SYMBOLS_SIDE_ISI + isiPlotXlimPadding]);
    ylim_rc = ylim(ax3);
    text(ax3, -NUM_SYMBOLS_SIDE_ISI, ylim_rc(2) * isiTextYPosFactor, ...
         isiTextContent, 'HorizontalAlignment', isiTextHAlign, ...
         'VerticalAlignment', isiTextVAlign, 'BackgroundColor', isiTextBGColor, ...
         'EdgeColor', isiTextEdgeColor, 'FontSize', isiTextFontSize);
    hold(ax3, 'off');
end