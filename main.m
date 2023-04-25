clear
%{
    Instructions
Choose no less than 2000 data collection length

    1001:
Fs = 5 MHz
TimeEnd = 20.736
    1002:
Fs = 5 MHz
TimeEnd = 
    1003:
Fs = 10 MHz
TimeEnd = 

%}
    %Define
experimentNo = '1001'; %Specify which experiment to analize
ASCIIOutPut = importdata(append('Data\EXP', experimentNo, '.txt'));
ASCIIWaveforms = append('Data\EXP', experimentNo);
ApplyHAF = false;
    %Hardware calibrations
PT = 20*10^-6; %Pre-trigger
PDT = 35; %Peak Definition Time
HDT = 150; %Hit Definition Time
HLT = 300; %Hit Lockout Time
Fs = 5*10^6; %Sample frequency (Hz)

    %Software parameters
Total = length(ASCIIOutPut.data)/Fs;
TimeEnd = 50; %Experiment cutoff time [s]
SampleNumber = 1; %Matrix crack no. to sample !OBS ERROR IF > TOTAL!
HAFfilter = 0; %Default: -1500

    %Matrixcrack definition
MCminFreq = 75*10^3; %[Hz]
MCmaxFreq = 180*10^3; %[Hz]
MCminAmp = 60; %[dB]
MCmaxAmp = 99; %[dB]
MCminEner = 0; %[kJ]
MCmaxEner = 10^12; %[aJ]
MCminDur = 1500; %[μs]

    %Debonding definition
DBminFreq = 240*10^3; %[Hz]
DBmaxFreq = 310*10^3; %[Hz]
DBminAmp = 45; %[dB]
DBmaxAmp = 65; %[dB]
DBminEner = 0; %[kJ]
DBmaxEner = 10^12; %[aJ]
DBminDur = 0; %[μs]

    %Import waveform
filePattern = fullfile(ASCIIWaveforms, '*.txt');
TheFiles = dir(filePattern);

    %Sort the files
C = transpose({TheFiles.name});
% Extract indices from each string using regular expressions
filesindices = cellfun(@(x) regexp(x, '(?<=_)\d+(?=_\d+\.txt)',...
    'match'), C, 'UniformOutput', false);
% Convert indices from cell array of strings to numeric array
filesindices = cellfun(@(x) str2double(x), filesindices, ...
    'UniformOutput', false);
% Sort the indices
[sorted_indices, order] = sort([filesindices{:}]);
% Reorder the original cell array using the sorted indices
TheFiles = TheFiles(order);

    %Find index of last sample given TimeEnd
FileTimes = str2double(regexp({TheFiles.name},...
    '(?<=_)\d+(?=\.txt)', 'match','once'))/1000000;
[~, HighestIndex] = max(FileTimes(FileTimes <= TimeEnd));

    %Import variables
DurationIndex = find(contains(ASCIIOutPut.textdata, 'DURATION'));
ImpDurList = transpose(ASCIIOutPut.data(:,DurationIndex));
ImpDurList = ImpDurList(1:HighestIndex);
AmplitudeIndex = find(contains(ASCIIOutPut.textdata, 'AMP'));
ImpAmpList = transpose(ASCIIOutPut.data(:,AmplitudeIndex));
ImpAmpList = ImpAmpList(1:HighestIndex);
EnergyIndex = find(contains(ASCIIOutPut.textdata, 'ENER'));
ImpEnerList = transpose(ASCIIOutPut.data(:,EnergyIndex));
ImpEnerList = ImpEnerList(1:HighestIndex);
PARA1Index = find(contains(ASCIIOutPut.textdata, 'PARA1'));
ImpPARA1 = transpose(ASCIIOutPut.data(:,PARA1Index));
ImpPARA1 = ImpPARA1(1:HighestIndex);

AMPList = [];
PFreqList = [];
HitTimeList = [];
HitIndexList = [];
Matrixcracks = [];
Debondings = [];
Delaminations = [];
MCc = 0; %MC counter
DBc = 0; %DB counter
FFTMat = [];
PowerMat = [];
SortEnerList = [];
StackEner = zeros(1,10^6);
SumFFT = zeros(1,HighestIndex);
BinFFT = zeros(1,HighestIndex);
HighAmpFilterHits = 0;
K = 0;
for k = 1 : HighestIndex
    K = K + 1;
        %Find file
    baseFileName = TheFiles(k).name;
    %baseFileName = FilesSorted(k);
    fullFileName = fullfile(TheFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', baseFileName);
        %Extract values
    Signal = load(fullFileName);
    N = length(Signal);  %number of samples    
    L = N/Fs; %total time of waveform (s)
        %Filter signal to duration
    if N > (PT+ImpDurList(k)*10^-6)*Fs %Check if duration > collection time
        FiltSignal = Signal(round(PT*Fs):round((PT+ImpDurList(k)*10^-6)*Fs));
    else
        FiltSignal = Signal(round(PT*Fs):N);
    end
    Nf = length(FiltSignal); %number of filtered samples
    Lf = Nf/Fs; %Total time of filtered waveform (s)
    dBpreamp = 0;
        %Operations
    HitTime = str2double(regexp(fullFileName,...
        '(?<=_)\d+(?=\.txt)', 'match','once'))/1000000; %Find time of hit
    HitIndex = str2double(regexp(fullFileName,...
        '(?<=_)\d+(?=_[^_]*\.txt)', 'match','once')); %Find index of hit
    
    if k == 1 %Define a spread vector/matrix
        Resolution = 20; %Recommended 1 - 20 Default: 10
        TimeTable = zeros(N,TimeEnd*Resolution);
        TimeVector = zeros(1,TimeEnd*Resolution);
        EnerList = TimeVector;
        SpreadEner = TimeVector;
        SpreadEner1 = TimeVector;
        SpreadEner2 = TimeVector;
        SpreadEner3 = TimeVector;
        DerivEner = TimeVector;
        DerivEner1 = TimeVector;
        DerivEner2 = TimeVector;
        DerivEner3 = TimeVector;
    end
    TimeIndex = round(HitTime*Resolution);

    FFTf = fft(FiltSignal); %Fast Fourier Transform
    FFT = fft(Signal);
    FFTMat(:,end+1) = abs(FFT(1:round(N/10))); %N for 10MHz limit
    for p = 1 : N
        SumFFT(k) = SumFFT(k)+log10(abs(FFT(p)));
    end
    if SumFFT(k) >= HAFfilter || ApplyHAF == false %Apply High Amplitude Filter
        HighAmpFilterHits = HighAmpFilterHits + 1;
        BinFFT(k) = 250;
        fVals = (0:Nf-1)/Lf;
        power = 20*log10(2*(abs(FFTf(1:round(Nf/2+1))))/Nf)+dBpreamp;
        [maxValue,indexMax] = max(abs(FFTf));
        PFreq = indexMax/Lf; %Peak frequency (Hz)
        AMP = max(power); %Amplitude
        PowerMat(:,end+1) = 20*log10(2*(abs(FFT(1:round(N/2+1))))/N);

        TimeTable(:,TimeIndex) = abs(FFT);

        SpreadEner(TimeIndex:end) = SpreadEner(TimeIndex:end)...
                                    +ImpEnerList(k);

        if 0 <= PFreq && PFreq <= 200*10^3 %Frequency band 1
            SpreadEner1(TimeIndex:end) = SpreadEner1(TimeIndex:end)...
                                        +ImpEnerList(k);
        end
        if 200*10^3 <= PFreq && PFreq <= 400*10^3 %Frequency band 2
            SpreadEner2(TimeIndex:end) = SpreadEner2(TimeIndex:end)...
                                        +ImpEnerList(k);
        end
        if 400*10^3 <= PFreq %Frequency band 3
            SpreadEner3(TimeIndex:end) = SpreadEner3(TimeIndex:end)...
                                        +ImpEnerList(k);
        end
        %disp(PFreq);
        StackEner(round(PFreq/1000)) = StackEner(round(PFreq/1000))...
                                     + ImpEnerList(k);

        PFreqList(K) = round(PFreq,1);
        AMPList(K) = round(AMP,1);
        HitTimeList(K) = HitTime;
        HitIndexList(K) = HitIndex;
        HAFImpAmpList(K) = ImpAmpList(k);
        HAFImpDurList(K) = ImpDurList(k);
        HAFImpEnerList(K) = ImpEnerList(k);
        HAFImpPARA1(K) = ImpPARA1(k);
            %Define matrix crack
        if MCminFreq <= PFreq && PFreq <= MCmaxFreq
            if MCminAmp <= ImpAmpList(k) && ImpAmpList(k) <= MCmaxAmp
                if MCminEner <= ImpEnerList(k) && ImpEnerList(k) <= MCmaxEner
                    if MCminDur <= ImpDurList(k)
                        MCc = MCc + 1;
                        if MCc == 1
                            Matrixcracks = table(); %Make table first time
                        end
                        if MCc == SampleNumber %Take a sample
                            SaFiltSignal = FiltSignal;
                            SaSignal = Signal;
                            SaDur = ImpDurList(k);
                            SaNf = Nf;
                            SaHitTime = HitTime;
                            SaHitIndex = HitIndex;
                            SaFFT = FFTf;
                            Sapower = power;
                            SafVals = fVals;
                            SaVals = (0:N-1)/L;
                            SaPFreq = PFreq;
                        end
                        Duration = ImpDurList(k);
                        Energy = ImpEnerList(k);
                        Amplitude = ImpAmpList(k);
                        PeakFrequency = PFreq/1000;
                        Load = ImpPARA1(k);
                        Matrixcracks(MCc,:) = table(HitIndex,HitTime,...
                            PeakFrequency,Amplitude,Duration,Energy,...
                            Load);
                    end
                end
            end
        end

            %Define debonding
        if DBminFreq <= PFreq && PFreq <= DBmaxFreq
            if DBminAmp <= ImpAmpList(k) && ImpAmpList(k) <= DBmaxAmp
                if DBminEner <= ImpEnerList(k) && ImpEnerList(k) <= DBmaxEner
                    if DBminDur <= ImpDurList(k)
                        DBc = DBc + 1;
                        if DBc == 1
                            Debondings = table(); %Make table first time
                        end
                        Duration = ImpDurList(k);
                        Energy = ImpEnerList(k);
                        Amplitude = ImpAmpList(k);
                        PeakFrequency = PFreq/1000;
                        Load = ImpPARA1;
                        Debondings(DBc,:) = table(HitIndex,HitTime,...
                            PeakFrequency,Amplitude,Duration,Energy,...
                            Load);
                    end
                end
            end
        end
    else
        K = K-1;
    end
end
%%
if ApplyHAF == true
    disp("HAF enabled");
    disp("Hits reduced:" + (HighestIndex - HighAmpFilterHits));
    disp("Hits remaining:" + HighAmpFilterHits);
else
    disp("HAF disabled");
end
for k = 1 : length(SpreadEner)-1
    DerivEner(k) = (SpreadEner(k+1)-SpreadEner(k))/(k+1);
    DerivEner1(k) = (SpreadEner1(k+1)-SpreadEner1(k))/(k+1);
    DerivEner2(k) = (SpreadEner2(k+1)-SpreadEner2(k))/(k+1);
    DerivEner3(k) = (SpreadEner3(k+1)-SpreadEner3(k))/(k+1);
end

close all
if istable(Debondings)
    disp("Debondings found:"...
        + num2str(length(Debondings.HitIndex)));
else
    disp("No debondings found");
end

if istable(Matrixcracks) %Sample waveforms
    disp("Matrix cracks found:"...
        + num2str(length(Matrixcracks.HitIndex)));
    figure;
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

figure('Position',[100 100 1000 500]); 
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

figure;
x = linspace(1,100,50);
tiledlayout(2,3);

nexttile %Amplitude-Frequency
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

nexttile %Hitcounter
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

nexttile %Total Accumulated Acoustic Energy
hold on
time2 = (0:length(TimeVector)-1);
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
legend("Total", "0-200kHz", "200-400kHz", ">400kHz",...
    'location','south outside');
hold off

nexttile %Energy-time derivative
hold on
plot(time2/Resolution, DerivEner);
plot(time2/Resolution, DerivEner1,'--');
plot(time2/Resolution, DerivEner2,'--');
plot(time2/Resolution, DerivEner3,'--');
xlabel('Time [s]');
ylabel('Energy derivative');
if istable(Matrixcracks) %If true = there are matrixcracks
    for i = 1 : length(Matrixcracks.HitTime)
        xli = xline(Matrixcracks.HitTime(i));
        xli.Color = [1 0 0];
    end
end
hold off

nexttile %Frequency vs Time vs Amplitude
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

%{
nexttile %Amplitude vs Time
hold on
scatter(HitTimeList, ImpAmpList, 60, '.');
refl = refline(0,MCminAmp); %minimum amp
refl.Color = 'r';
refl = refline(0,MCmaxAmp); %maximum amp
refl.Color = 'r';
title('Amplitude vs Time');
%xlim([0 TimeEnd]);
ylim([floor(min(ImpAmpList))*0.9 ceil(max(ImpAmpList))*1.1]);
xlabel('Time [s]');
ylabel("Amplitude [dB]");
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.HitTime, Matrixcracks.Amplitude, 'o');
end
hold off
%}

nexttile %Energy-Frequency Spectrum
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

%Primary 2
figure;
x = linspace(1,100,50);
tiledlayout(2,3);

nexttile %Duration vs Time
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

nexttile %Energy vs Time
hold on
scatter(HitTimeList, HAFImpEnerList, 60, '.');
title('Energy vs time');
%xlim([0 TimeEnd]);
ylim([0 3*10^5]);
xlabel('Time [s]');
ylabel("Energy [aJ]");
refl = refline(0,MCminEner); %minimum duration
refl.Color = 'r';
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.HitTime, Matrixcracks.Energy, 'o');
end
hold off

nexttile %Load

hold on
scatter(HAFImpPARA1, HAFImpAmpList, 60, '.');
%refl = refline(0,MCminAmp); %minimum amp
%refl.Color = 'r';
%refl = refline(0,MCmaxAmp); %maximum amp
%refl.Color = 'r';
title('Amplitude vs Load');
%xlim([0 TimeEnd]);
ylim([floor(min(HAFImpAmpList))*0.9 ceil(max(HAFImpAmpList))*1.1]);
xlabel('Load [ ]');
ylabel("Amplitude [dB]");
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.Load, Matrixcracks.Amplitude, 'o');
end
hold off
%{
nexttile %Frequency vs Time vs Energy
hold on
PeakFrequencyList = PFreqList/1000;
tbl = table(HitTimeList,PeakFrequencyList,ImpEnerList);
sca = scatter(HitTimeList,PeakFrequencyList,60,ImpEnerList,'.');
sca.CData = tbl.ImpEnerList;
cb = colorbar;
refl = refline(0,MCminFreq/1000); %minimum frequency
refl.Color = 'r';
refl = refline(0,MCmaxFreq/1000); %maximum frequency
refl.Color = 'r';
title(cb,'dB')
title('Frequency vs Time vs Energy');
%xlim([0 TimeEnd]);
ylim([0 1000]);
xlabel('Time [s]');
ylabel("Peak Frequency [kHz]");
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.HitTime, Matrixcracks.PeakFrequency, 'o');
end
hold off
%}

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