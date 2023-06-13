function data = preprocess_fn(filename, savepath, tfdana)
    tfdana = 1; % 1 provide time-frequency estimate, otherwise 0

    % LOAD DATA
    Signals = load(filename);
    IDch = {'1Y' '1X' '2Y'};

    sign_IDCh = [-1 1, -1];
    tdata = double([Signals.Data1_AI_2 Signals.Data1_AI_3, Signals.Data1_AI_4]); % acc in [m/s2]
    fs = Signals.Sample_rate; % [Hz]

    % desampling of the signals -----------------------------------------------
    RR = 1; % desampling factor
    for ii = 1:size(tdata, 2)
        tdataDES(:, ii) = decimate(tdata(:, ii), RR);
    end
    fs = fs / RR;
    tdata = tdataDES;
    tdataDES = [];
    % -------------------------------------------------------------------------

    IDRetCh = {'1Y' '1X' '2Y'};
    RetCh = 1:size(tdata, 2); % retained sensor channels --> analysis was performed exploring the time and frequency series of the data
    tdata = tdata(:, RetCh);
    IDch = IDch(RetCh);
    sign_IDCh = sign_IDCh(RetCh);

    % PROCESSING DATA
    Output = sign_IDCh .* tdata;
    Output = Output - nanmean(Output);

    % -------------------------------------------------------------------------

    % Calculating spectral quantities
    ts = 1 / fs;
    time = ts * (0:size(Output, 1) - 1)';
    tend = time(end);
    fres = 1 / tend;
    freq = -fs / 2:fres:fs / 2; %(0:length(time)-1)'*fres-fs/2;

    fcutmin = 0.5; % [Hz]
    fcutmax = 15; % [Hz]
    [BB, AA] = butter(4, [fcutmin fcutmax] / (fs / 2));
    Output = filtfilt(BB, AA, Output);

    OUTPUT = fftshift(fft(Output, [], 1) * ts, 1);
    OUTPUTa = abs(OUTPUT);
    FourSpec = OUTPUTa .^ 2;

    [Pwelch, Freq] = pwelch(Output, [], [], length(Output), fs);

    %TIME-FREQUENCY ANALYSIS
    mag_fac = 1; % magnification factor for time-frequency plot

    if tfdana == 1 %1 solo quando ti serve controllare o vedere ?
        Powe = 15; % raw wstimate of freqres & winlengthS
        freqres = 1 / (2^Powe * ts); % frequency resolution [Hz] of the TFD
        winlengthS = 1 / freqres; % time length of the window [s] for the TFD
        winlengthN = round(winlengthS / ts); % discrete length of the window [samples] for the TFD
        win = hamming(winlengthN, 'periodic'); % window (periodic because the spectral analysis)
        overlap = 99 / 100; % win overlap

        OutputTFDcspecNormSum = 0;
        for ch = 1:size(Output, 2)
            if ch == 1
                [~, FRE_TFD, TIM_TFD, OutputTFDc{ch}] = spectrogram(Output(:, ch), win, round(overlap * (length(win) - 1)), 2 * round(fs / freqres), fs, 'psd');
            else
                [~, ~, ~, OutputTFDc{ch}] = spectrogram(Output(:, ch), win, round(overlap * (length(win) - 1)), 2 * round(fs / freqres), fs, 'psd');
            end
            OutputTFDcspec = OutputTFDc{ch} .^ mag_fac;
            OutputTFDcspecNorm = OutputTFDcspec ./ max(abs(OutputTFDcspec), [], 'all');
            OutputTFDcspecNormSum = OutputTFDcspecNormSum + OutputTFDcspecNorm;
        end
        OutputTFDcspecNormMean = (OutputTFDcspecNormSum / size(Output, 2));
        OutputTFDcspecNormMeanFREQ = nanmean(OutputTFDcspecNormMean, 2);
    end

    data = [time Output];
    save(savepath, 'data');
end

