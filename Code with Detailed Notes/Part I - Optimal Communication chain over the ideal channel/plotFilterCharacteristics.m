function hFig = plotFilterCharacteristics(h_rrc, Beta, Fs, OSF)
    %PLOTFILTERCHARACTERISTICS Visualizes RRC and combined RC filter properties.
    %   HFIG = PLOTFILTERCHARACTERISTICS(H_RRC, BETA, FS, OSF) creates a figure
    %   with three subplots analyzing the provided Root Raised Cosine (RRC)
    %   filter impulse response (H_RRC):
    %   1. Top Plot: RRC filter impulse response in the time domain.
    %   2. Middle Plot: RRC filter magnitude frequency response (dB scale).
    %   3. Bottom Plot: Combined Raised Cosine (RC) pulse response in the time
    %      domain (obtained by convolving H_RRC with itself), normalized by
    %      symbol period. Includes markers at symbol sampling instants (kT_s)
    %      to visually check for Inter-Symbol Interference (ISI).
    %
    %   Assumes H_RRC is a valid impulse response vector. Plotting appearance
    %   is controlled by internal parameters.
    %
    %   Inputs:
    %       h_rrc             - Vector containing the RRC filter coefficients
    %                           (time-domain impulse response).
    %       Beta              - The roll-off factor (0 <= Beta <= 1) used to
    %                           design the RRC filter. Used for annotation.
    %       Fs                - The sampling frequency (in Hz) at which the
    %                           filter H_RRC operates. Used for frequency axis scaling.
    %       OSF               - The oversampling factor (samples per symbol)
    %                           used in the system. Used for ISI check calculation.
    %
    %   Output:
    %       hFig              - Handle to the created figure object.

    % =====================================================================
    % == Plotting Parameters (Locally Defined) ==
    % =====================================================================
    % --- Colors ---
    filterColor = [0.4660, 0.6740, 0.1880]; % Greenish color for RRC filter plots
    linkColor   = [0.4940, 0.1840, 0.5560]; % Purple color for the combined RC pulse plot
    rxColor     = [0.8500, 0.3250, 0.0980]; % Reddish color for ISI check markers
    % --- Line/Marker Styles ---
    lineWidth = 1.5;                        % Line width for continuous plots
    baseMarkerSize = 6;                     % Base size for stem markers
    rrcImpulseStemStyle = 'filled';         % Style for RRC impulse response stem plot ('filled', 'open')
    rrcImpulseMarkerSize= 4;                % Marker size for RRC impulse stems
    isiStemMarkerShape  = 'o';              % Marker shape for ISI check stems ('o', 'x', etc.)
    isiStemLineStyle    = 'none';           % Line style for ISI stems ('-', 'none', etc.)
    % --- General Figure Settings ---
    figureName          = 'Filter Analysis (RRC & RC)'; % Figure window title
    figurePosition      = [100, 100, 800, 700]; % Position and size [left, bottom, width, height]
    useFigureNumberTitle= 'off';            % Hide 'Figure X:' prefix ('on' or 'off')
    % --- Calculation Parameters ---
    NFFT                = 2^12;             % Number of FFT points for frequency response calculation (>= NumTaps)
    NUM_SYMBOLS_SIDE_ISI= 5;                % Number of symbol periods to the left/right of center for ISI check plot
    % --- Unit Conversions ---
    timeUnitsFactor     = 1e6;              % Factor to convert seconds to microseconds (us)
    timeUnitsLabel      = '\mus';           % Label for microsecond units
    freqUnitsFactor     = 1e6;              % Factor to convert Hz to Megahertz (MHz)
    freqUnitsLabel      = 'MHz';            % Label for Megahertz units
    % --- Subplot 1: RRC Impulse Response ---
    subplot1TitleFormat = 'RRC Filter Impulse Response (β=%.2f, Taps=%d)'; % Title format
    subplot1XLabel      = sprintf('Time (%s)', timeUnitsLabel); % X-axis label
    subplot1YLabel      = 'Amplitude';          % Y-axis label
    % --- Subplot 2: RRC Frequency Response ---
    subplot2Title       = 'RRC Filter Frequency Response (Magnitude)'; % Title
    subplot2XLabel      = sprintf('Frequency (%s)', freqUnitsLabel); % X-axis label
    subplot2YLabel      = 'Magnitude (dB)';     % Y-axis label
    freqResponseMinDb   = -80;              % Lower limit for y-axis (dB) in freq response plot
    freqResponseMaxDb   = 10;               % Upper limit for y-axis (dB) in freq response plot
    % --- Subplot 3: Combined RC Response & ISI Check ---
    subplot3Title       = 'Combined Raised Cosine Pulse Response & ISI Check'; % Title
    subplot3XLabel      = 'Time (Normalized to Symbol Periods, t/T_{symb})'; % X-axis label (normalized time)
    subplot3YLabel      = 'Amplitude';          % Y-axis label
    isiPlotXlimPadding  = 0.5;              % Padding for x-axis limits in ISI plot (in symbol periods)
    isiLegendLocation   = 'best';           % Location for the legend in the ISI plot
    isiTextContent      = ' Samples at k*T_{symb} '; % Text annotation for ISI markers
    isiTextYPosFactor   = 0.9;              % Vertical position factor for ISI text (relative to ylim)
    isiTextHAlign       = 'left';           % Horizontal alignment for ISI text
    isiTextVAlign       = 'top';            % Vertical alignment for ISI text
    isiTextBGColor      = 'w';              % Background color for ISI text box
    isiTextEdgeColor    = 'k';              % Edge color for ISI text box
    isiTextFontSize     = 8;                % Font size for ISI text
    % =====================================================================

    % =====================================================================
    % == Derived Parameters & Constants ==
    % =====================================================================
    % Get the number of taps from the input impulse response vector.
    NumFilterTaps = length(h_rrc);
    % Calculate the sample period (duration of one sample).
    Ts = 1/Fs; % Sample period in seconds
    % Calculate the symbol period based on sample period and OSF.
    Tsymb = Ts * OSF; % Symbol period in seconds

    % =====================================================================
    % == Plotting Code ==
    % =====================================================================
    % Create the main figure window with specified properties.
    hFig = figure('Name', figureName, ...
                  'NumberTitle', useFigureNumberTitle, ...
                  'Position', figurePosition);

    % --- Subplot 1: RRC Impulse Response (Time Domain) ---
    ax1 = subplot(3, 1, 1, 'Parent', hFig); % Create the first subplot
    % Calculate the time axis for the filter taps, centered around t=0.
    time_axis_filter_sec = (-(NumFilterTaps-1)/2 : (NumFilterTaps-1)/2) * Ts;
    % Convert time axis to desired units (e.g., microseconds).
    time_axis_filter_units = time_axis_filter_sec * timeUnitsFactor;
    % Plot the RRC impulse response using a stem plot.
    stem(ax1, time_axis_filter_units, h_rrc, ...
        rrcImpulseStemStyle, ...            % Use specified stem style
        'Color', filterColor, ...           % Use specified color
        'MarkerSize', rrcImpulseMarkerSize); % Use specified marker size
    grid(ax1, 'on'); box(ax1, 'on');       % Enable grid and box
    xlabel(ax1, subplot1XLabel);           % Set x-axis label
    ylabel(ax1, subplot1YLabel);           % Set y-axis label
    title(ax1, sprintf(subplot1TitleFormat, Beta, NumFilterTaps)); % Set title with parameters
    xlim(ax1, time_axis_filter_units([1 end])); % Adjust x-limits to fit the data

    % --- Subplot 2: RRC Frequency Response (Magnitude) ---
    ax2 = subplot(3, 1, 2, 'Parent', hFig); % Create the second subplot
    % Calculate the frequency response using freqz.
    % 'whole' computes the response over [0, 2*pi) or [0, Fs).
    [H_full, ~] = freqz(h_rrc, 1, NFFT, 'whole', Fs);
    % Shift the frequency response so that DC (0 Hz) is in the center.
    H_shifted = fftshift(H_full);
    % Create the corresponding frequency axis, centered at 0 Hz.
    F_shifted_hz = (-NFFT/2 : NFFT/2 - 1) * (Fs / NFFT);
    % Convert frequency axis to desired units (e.g., MHz).
    F_shifted_units = F_shifted_hz / freqUnitsFactor;
    % Plot the magnitude response in dB. Add 'eps' for numerical stability (avoid log10(0)).
    plot(ax2, F_shifted_units, 20*log10(abs(H_shifted) + eps), ...
        'Color', filterColor, ...           % Use specified color
        'LineWidth', lineWidth);            % Use specified line width
    grid(ax2, 'on'); box(ax2, 'on');       % Enable grid and box
    xlabel(ax2, subplot2XLabel);           % Set x-axis label
    ylabel(ax2, subplot2YLabel);           % Set y-axis label
    title(ax2, subplot2Title);             % Set title
    xlim(ax2, [-Fs/2, Fs/2] / freqUnitsFactor); % Set x-limits to [-Fs/2, Fs/2] in target units
    ylim(ax2, [freqResponseMinDb, freqResponseMaxDb]); % Set y-limits for dB scale

    % --- Subplot 3: Combined RC Response & ISI Check ---
    ax3 = subplot(3, 1, 3, 'Parent', hFig); % Create the third subplot
    % Calculate the impulse response of the combined Raised Cosine (RC) filter.
    % This represents the end-to-end pulse shape seen at the receiver after
    % matched filtering, assuming the transmitter used the same RRC filter.
    % It's obtained by convolving the RRC impulse response with itself.
    h_rc = conv(h_rrc, h_rrc);
    % Determine the number of taps in the resulting RC filter response.
    NumRCTaps = length(h_rc);
    % Find the index corresponding to the center (peak) of the RC pulse.
    rc_center_index = floor(NumRCTaps / 2) + 1;
    % Calculate the time axis for the RC pulse, centered around t=0.
    time_axis_rc_sec = (-(NumRCTaps-1)/2 : (NumRCTaps-1)/2) * Ts;
    % Normalize the time axis by the symbol period (Tsymb) for ISI visualization.
    time_axis_rc_norm = time_axis_rc_sec / Tsymb;
    % Plot the continuous-time RC pulse shape.
    plot(ax3, time_axis_rc_norm, h_rc, ...
        'Color', linkColor, ...             % Use specified color for RC pulse
        'LineWidth', lineWidth);            % Use specified line width
    hold(ax3, 'on'); % Hold plot for adding ISI markers

    % --- ISI Check: Mark samples at integer multiples of Tsymb ---
    % Define the range of symbol intervals around the center (k=0) to check.
    symbol_indices_offset = -NUM_SYMBOLS_SIDE_ISI:1:NUM_SYMBOLS_SIDE_ISI; % e.g., -5, -4, ..., 0, ..., 4, 5
    % Calculate the sample indices corresponding to these symbol instants (k * Tsymb).
    % We need to offset from the center index of the h_rc vector. Each symbol period Tsymb corresponds to OSF samples.
    symbol_sample_indices = rc_center_index + symbol_indices_offset * OSF;
    % Ensure indices are valid (within the bounds of h_rc).
    valid_indices_mask = (symbol_sample_indices >= 1) & (symbol_sample_indices <= NumRCTaps);
    symbol_sample_indices = symbol_sample_indices(valid_indices_mask);
    symbol_indices_offset = symbol_indices_offset(valid_indices_mask); % Keep corresponding offsets
    % Get the normalized time values at these symbol sampling instants.
    symbol_times_norm = time_axis_rc_norm(symbol_sample_indices);
    % Get the amplitude of the RC pulse at these symbol sampling instants.
    symbol_values = h_rc(symbol_sample_indices);
    % Plot stems at these points to visualize ISI. Ideally, values should be non-zero
    % only at k=0 and zero (or very small) for k != 0.
    stem(ax3, symbol_times_norm, symbol_values, ...
         'Color', rxColor, ...              % Use specified color for markers
         'Marker', isiStemMarkerShape, ...  % Use specified marker shape
         'MarkerFaceColor', rxColor, ...    % Fill marker with the same color
         'MarkerSize', baseMarkerSize, ...  % Use specified marker size
         'LineStyle', isiStemLineStyle);    % Use specified line style for stems
    % Add a legend to clarify the plot elements.
    legend(ax3, 'RC Pulse Shape', sprintf('Samples at k*T_{symb} (OSF=%d)', OSF), ...
           'Location', isiLegendLocation);
    grid(ax3, 'on'); box(ax3, 'on');       % Enable grid and box
    xlabel(ax3, subplot3XLabel);           % Set x-axis label (normalized time)
    ylabel(ax3, subplot3YLabel);           % Set y-axis label
    title(ax3, subplot3Title);             % Set title
    % Set x-axis limits to focus on the specified number of symbol periods around the center.
    xlim(ax3, [-NUM_SYMBOLS_SIDE_ISI - isiPlotXlimPadding, NUM_SYMBOLS_SIDE_ISI + isiPlotXlimPadding]);
    % Add text annotation near the ISI markers for clarity.
    ylim_rc = ylim(ax3); % Get current y-axis limits to position text
    text(ax3, -NUM_SYMBOLS_SIDE_ISI, ylim_rc(2) * isiTextYPosFactor, ... % Position text near top-left of ISI region
         isiTextContent, ...
         'HorizontalAlignment', isiTextHAlign, ...
         'VerticalAlignment', isiTextVAlign, ...
         'BackgroundColor', isiTextBGColor, ...
         'EdgeColor', isiTextEdgeColor, ...
         'FontSize', isiTextFontSize);
    hold(ax3, 'off'); % Release the plot hold
end