if istable(Matrixcracks) %Sample waveforms
    disp("Matrix cracks found:"...
        + num2str(length(Matrixcracks.HitIndex)));
    figure('name', 'Sample Waveform')
    tiledlayout(2,2);

    nexttile %Waveform
    plot(SaVals*L/Fs*10^6, SaSignal);
    title('Sample Waveform');
    xlabel('Time [μs]');
    ylabel('Voltage [V]');
    refline(0,max(SaSignal)*2); %top threshold
    refline(0,-max(SaSignal)*2); %bot threshold
    xline(PT*10^6); %Pretrigger line
    xline(PT*10^6+SaDur); %Filter line
    xlim([0 max(SaVals*L/Fs*10^6)]);
    %ylim([-ceil(max(FiltSignal)), ceil(max(FiltSignal))]);

    nexttile %Power spectrum
    plot(SafVals(1:round(SaNf/2+1)), Sapower);
    xlim([0 10^6]);
    title('Sample Power spectrum');
    xlabel('Frequency [Hz]');
    ylabel('Amplitude [dB]');
    example = 'example_string_123';

    nexttile %FFT
    plot(SafVals(1:round(SaNf/2+1)), abs(SaFFT(1:round(SaNf/2+1))));
    title('Sample FFT');
    xlim([0 10^6]);
    xlabel('Frequency [Hz]');
    ylabel('Voltage [V]');
else
        disp("No matrix cracks found");
end

figure('name', 'Spectrogram', 'Position',[100 100 1000 500]); 
tiledlayout(2,1);
nexttile %Spectogram
image([0 HighestIndex], [0 1*10^6], FFTMat);
title('Spectrogram');
cb = colorbar;
title(cb,'        Intensity')
xlabel('Waveform no.');
ylabel('Frequency [kHz]');
annotation('textbox', [0.8 0.87 0.8 0.1], ...
    'String', append('Total hits:', num2str(HighestIndex)), ...
    'Color', [1 0 0], ...
    'FontWeight', 'bold', ...
    'EdgeColor', 'none')

nexttile %High amplitude filter
image([0 HighestIndex], [0 0], BinFFT);
title(append('HAF (',num2str(HAFfilter),')'));
xlabel('Waveform no.');
annotation('textbox', [0.79 0.393 0.8 0.1], ...
    'String', append('Filtered hits:', num2str(K)), ...
    'Color', [1 0 0], ...
    'FontWeight', 'bold', ...
    'EdgeColor', 'none')

figure('name', 'Amplitude-Frequency')
hold on
plot(PFreqList./1000, HAFImpAmpList, 'x');
xlim([floor(min(PFreqList./1000)/10)*10 ...
    ceil(max(PFreqList./1000)/10)*10]); %nearest 10-number of max sample
ylim([floor(min(HAFImpAmpList)/10)*10 ...
    ceil(max(HAFImpAmpList)/10)*10]);
title('Amplitude-Frequency');
xlabel('Peak frequency [kHz]');
ylabel('Amplitude [dB]');
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.PeakFrequency, Matrixcracks.Amplitude, 'o');
    legend(plus("Hits: ",num2str(K)),...
        plus("Matrix cracks: ",...
        num2str(length(Matrixcracks.HitIndex))),...
        'location','south outside',"AutoUpdate","off");
else
    legend(plus("Hits: ",num2str(K)),...
        plus("Matrix cracks: ",...
        num2str(length(Matrixcracks))),'location','south outside',"AutoUpdate","off");
end
patch([MCminFreq/1000 MCminFreq/1000 MCmaxFreq/1000 MCmaxFreq/1000],...
    [MCminAmp MCmaxAmp MCmaxAmp MCminAmp],'r','FaceAlpha',0,...
    'EdgeColor','r');
hold off

figure('name', 'Hitcounter')
hold on
%time = (0:ceil(max(HitTimeList)));
time = (0:TimeEnd);
SpreadHits = zeros(1,length(time));
for i = 1 : length(HitTimeList)
    SpreadHits(ceil(HitTimeList(i)):end) = ...
        SpreadHits(ceil(HitTimeList(i)):end)+1;
end
plot(time, SpreadHits);
title('Hits over time');
xlim([0 max(HitTimeList)]);
ylim([0 max(SpreadHits)+2]);
xlabel('Time [s]');
ylabel('Hits');
xline(PullStop);
DBSpreadHits = zeros(1,length(time)); %Debondings
if istable(Debondings) %If true = there are matrixcracks
    for i = 1 : length(Debondings.HitIndex)
        DBSpreadHits(ceil(Debondings.HitTime(i)):end) = ...
            DBSpreadHits(ceil(Debondings.HitTime(i)):end)...
            + 1;
    end
end
if istable(Debondings) %If true = there are debondings
    plot(time, DBSpreadHits,'g');
end
MCSpreadHits = zeros(1,length(time)); %Matrix cracks
if istable(Matrixcracks) %If true = there are matrixcracks
    for i = 1 : length(Matrixcracks.HitIndex)
        MCSpreadHits(ceil(Matrixcracks.HitTime(i)):end) = ...
            MCSpreadHits(ceil(Matrixcracks.HitTime(i)):end)...
            + 1;
    end
end
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(time, MCSpreadHits,'r');
end
legend("Hits","Debondings","Matrix cracks",'location','south outside');
hold off

figure('name', 'Total Accumulated Acoustic Energy')
hold on
plot(time2/Resolution, SpreadEner);
% FOR POLYNOMIAL % 
% p = polyfit(time2/Resolution,SpreadEner,20);
% pol = polyval(p, time2/Resolution);
% plot(time2/Resolution, pol,'--');
title('Accumulated Acoustic Energy');
xlim([0 TimeEnd]);
%ylim([0 max(SpreadEner)+2]);
xlabel('Time [s]');
ylabel('Energy');
plot(time2/Resolution, SpreadEner1,'--');
plot(time2/Resolution, SpreadEner2,'--');
plot(time2/Resolution, SpreadEner3,'--');
xline(PullStop);
legend("Total", "0-200kHz", "200-400kHz", ">400kHz",...
    'location','south outside');
hold off

figure('name', 'Energy-time derivative')
hold on
plot(time2/Resolution, DerivEner);
plot(time2/Resolution, DerivEner1,'--');
plot(time2/Resolution, DerivEner2,'--');
plot(time2/Resolution, DerivEner3,'--');
title('Energy Derivative vs Time');
xlabel('Time [s]');
ylabel('Energy derivative');
if istable(Matrixcracks) %If true = there are matrixcracks
    for i = 1 : length(Matrixcracks.HitTime)
        xli = xline(Matrixcracks.HitTime(i));
        xli.Color = [1 0 0];
    end
end
hold off

figure('name', 'Frequency vs Time vs Amplitude')
hold on
PeakFrequencyList = PFreqList/1000;
tbl = table(HitTimeList,PeakFrequencyList,HAFImpAmpList);
sca = scatter(HitTimeList,PeakFrequencyList,60,HAFImpAmpList,'.');
sca.CData = tbl.HAFImpAmpList;
cb = colorbar;
refl = refline(0,MCminFreq/1000); %minimum frequency
refl.Color = 'r';
refl = refline(0,MCmaxFreq/1000); %maximum frequency
refl.Color = 'r';
title(cb,'dB')
title('Frequency vs Time vs Amplitude');
%xlim([0 TimeEnd]);
ylim([0 1000]);
xlabel('Time [s]');
ylabel("Peak Frequency [kHz]");
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.HitTime, Matrixcracks.PeakFrequency, 'o');
end
hold off

figure('name', 'Energy-Frequency Spectrum')
hold on
plot((1:length(StackEner)),StackEner);
xlim([0 1000]); 
ylim([0 max(StackEner)]);
title('Energy-Frequency Spectrum');
xlabel('Peak frequency [kHz]');
ylabel('Energy [aJ]');
xl1 = xline(MCminFreq/1000);
xl1.Color = 'r';
xl2 = xline(MCmaxFreq/1000);
xl2.Color = 'r';

if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.PeakFrequency, Matrixcracks.Energy, 'o');
end
hold off

figure('name', 'Duration vs Time')
hold on
scatter(HitTimeList, HAFImpDurList, 60, '.');
title('Duration');
%xlim([0 TimeEnd]);
xlabel('Time [s]');
ylabel("Duration [μs]");
refl = refline(0,MCminDur); %minimum duration
refl.Color = 'r';
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.HitTime, Matrixcracks.Duration, 'o');
end
hold off

figure('name', 'Energy vs Time')
hold on
scatter(HitTimeList, HAFImpEnerList, 60, '.');
title('Energy vs time');
ylim([0 3*10^5]);
xlabel('Time [s]');
ylabel("Energy [aJ]");
refl = refline(0,MCminEner); %minimum duration
refl.Color = 'r';
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.HitTime, Matrixcracks.Energy, 'o');
end
hold off

figure('name','Amplitude vs Load')
hold on
scatter(HAFImpPARA1, HAFImpAmpList, 60, '.');
title('Amplitude vs Load');
ylim([floor(min(HAFImpAmpList))*0.9 ceil(max(HAFImpAmpList))*1.1]);
xlabel('Load [ ]');
ylabel("Amplitude [dB]");
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.Load, Matrixcracks.Amplitude, 'o');
end
hold off

figure('name', 'Load')
hold on
plot(CSVDataOffs.Fun_Time, CSVDataOffs.Fun_Load)
hold off

if istable(Matrixcracks)
        %Compute statistics
    MCAvgFreq = round(mean(Matrixcracks.PeakFrequency),1);
    MCAvgAmp = round(mean(Matrixcracks.Amplitude),1);
    MCAvgDur = round(mean(Matrixcracks.Duration),1);
    MCAvgEner = round(mean(Matrixcracks.Energy),1);

    MCDevFreq = max(Matrixcracks.PeakFrequency) - round(MCAvgFreq,1);
    MCDevAmp = max(Matrixcracks.Amplitude) - round(MCAvgAmp,1);
    MCDevDur = max(Matrixcracks.Duration) - round(MCAvgDur,1);
    MCDevEner = max(Matrixcracks.Energy) - round(MCAvgEner,1);
    
    disp("MATRIX CRACK STATISTICS");
    disp("Avg peak frequency: " + num2str(MCAvgFreq) + "kHz ± " + num2str(MCDevFreq)); 
    disp("Avg amplitude: " + num2str(MCAvgAmp) + "dB ± " + num2str(MCDevAmp)); 
    disp("Avg duration: " + num2str(MCAvgDur) + "μs ± " + num2str(MCDevDur)); 
    disp("Avg energy: " + num2str(MCAvgEner) + "aJ ± " + num2str(MCDevEner)); 
end
if istable(Debondings)
    DBAvgFreq = round(mean(Debondings.PeakFrequency),1);
    DBAvgAmp = round(mean(Debondings.Amplitude),1);
    DBAvgDur = round(mean(Debondings.Duration),1);
    DBAvgEner = round(mean(Debondings.Energy),1);

    DBDevFreq = max(Debondings.PeakFrequency) - round(DBAvgFreq,1);
    DBDevAmp = max(Debondings.Amplitude) - round(DBAvgAmp,1);
    DBDevDur = max(Debondings.Duration) - round(DBAvgDur,1);
    DBDevEner = max(Debondings.Energy) - round(DBAvgEner,1);
    
    disp("DEBONDING STATISTICS");
    disp("Avg peak frequency: " + num2str(DBAvgFreq) + "kHz ± " + num2str(DBDevFreq)); 
    disp("Avg amplitude: " + num2str(DBAvgAmp) + "dB ± " + num2str(DBDevAmp)); 
    disp("Avg duration: " + num2str(DBAvgDur) + "μs ± " + num2str(DBDevDur)); 
    disp("Avg energy: " + num2str(DBAvgEner) + "aJ ± " + num2str(DBDevEner)); 
end
disp('Done.');