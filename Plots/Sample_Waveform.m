if istable(Matrixcracks) %Sample waveforms
    disp("Matrix cracks found:"...
        + num2str(length(Matrixcracks.HitIndex)));
    figure('name', 'Sample Waveform','Position',[60,60,1400,700]) %Waveform
    plot(SaVals*L/Fs*10^6, SaSignal);
    title('Sample Waveform');
    xlabel('Time [μs]');
    ylabel('Voltage [V]');
    refline(0,max(SaSignal)*2); %top threshold
    refline(0,-max(SaSignal)*2); %bot threshold
    xline(PT*10^6); %Pretrigger line
    xline(PT*10^6+SaDur); %Filter line
    xlim([0 max(SaVals*L/Fs*10^6)]);
    set(gca,'FontSize',14)
    %ylim([-ceil(max(FiltSignal)), ceil(max(FiltSignal))]);

    figure('name', 'Power spectrum','Position',[60,60,1400,700]) %Power spectrum
    plot(SafVals(1:round(SaNf/2+1)), Sapower);
    xlim([0 10^6]);
    title('Sample Power spectrum');
    xlabel('Frequency [Hz]');
    ylabel('Amplitude [dB]');
    example = 'example_string_123';
    set(gca,'FontSize',14)
    
    figure('name', 'FFT','Position',[60,60,1400,700]) %FFT
    area(SafVals(1:round(SaNf/2+1)), abs(SaFFT(1:round(SaNf/2+1))));
    title('Sample FFT');
    xlim([0 10^6]);
    xlabel('Frequency [Hz]');
    ylabel('Voltage [V]');
    legend('Total integral amplitude');
    set(gca,'FontSize',14)
else
        disp("No matrix cracks found");
end