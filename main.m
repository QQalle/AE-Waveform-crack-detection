close all
clear

    %Instructions
%Choose no less than 2000 data collection length

    %Define
experimentNo = '1002'; %Specify which experiment to analize
ASCIIOutPut = importdata(append('Data\EXP', experimentNo, '.txt'));
ASCIIWaveforms = append('Data\EXP', experimentNo);

    %Hardware calibrations
PT = 20*10^-6; %Pre-trigger
PDT = 35; %Peak Definition Time
HDT = 150; %Hit Definition Time
HLT = 300; %Hit Lockout Time
Fs = 5*10^6; %Sample frequency (Hz)

    %Software parameters
Total = length(ASCIIOutPut.data)/Fs;
TimeEnd = 80; %Experiment cutoff time [s]
SampleNumber = 2; %Matrix crack number in order of happening

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
%TimeIndex = find(contains(ASCIIOutPut.textdata, 'SSSSSSSS.mmmuuun'));
%ImpTimeList = transpose(ASCIIOutPut.data(:,TimeIndex));
%ImpTimeList = ImpTimeList(1:HighestIndex);
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
StackEner = [];
for k = 1 : HighestIndex
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
        Resolution = 2000;
        TimeTable = zeros(N,TimeEnd*Resolution);
        TimeVector = zeros(1,TimeEnd*Resolution);
        EnerList = TimeVector;
        SpreadEner = TimeVector;
        SpreadEner1 = TimeVector;
        SpreadEner2 = TimeVector;
        SpreadEner3 = TimeVector;
    end
    TimeIndex = round(HitTime*Resolution);
    
    FFTf = fft(FiltSignal); %Fast Fourier Transform
    FFT = fft(Signal);
    fVals = (0:Nf-1)/Lf;
    power = 20*log10(2*(abs(FFTf(1:round(Nf/2+1))))/Nf)+dBpreamp;
    [maxValue,indexMax] = max(abs(FFTf));
    PFreq = indexMax/Lf; %Peak frequency (Hz)
    AMP = max(power); %Amplitude
    FFTMat(:,end+1) = FFT;
    PowerMat(:,end+1) = 20*log10(2*(abs(FFT(1:round(N/2+1))))/N);
    %SpecMat(:,end+1) = spectrogram(Signal);
        %Add to list
    
    
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
    
    %StackEner(round(PFreqList(k))) = 
    
    PFreqList(k) = round(PFreq,1);
    AMPList(k) = round(AMP,1);
    HitTimeList(k) = HitTime;
    HitIndexList(k) = HitIndex;
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
end
%%
close all
if istable(Debondings)
    disp("Debondings found:"...
        + num2str(length(Debondings.HitIndex)));
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
    %xlim([0 Fs]);
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
    %xlim([0 10^6]);
    xlabel('Frequency [Hz]');
    ylabel('Voltage [V]');
else
        disp("No matrix cracks found");
end
%a = abs(SaFFT(1:Nf/2+1));
    %Primary figure
figure;
x = linspace(1,100,50);
tiledlayout(2,3);

nexttile %Amplitude-Frequency Spectrum
hold on
plot(PFreqList./1000, ImpAmpList, 'x');
xlim([floor(min(PFreqList./1000)/10)*10 ...
    ceil(max(PFreqList./1000)/10)*10]); %nearest 10-number of max sample
ylim([floor(min(ImpAmpList)/10)*10 ...
    ceil(max(ImpAmpList)/10)*10]);
title('Amplitude-Frequency Spectrum');
xlabel('Peak frequency [kHz]');
ylabel('Amplitude [dB]');
xl1 = xline(MCminFreq/1000);
xl1.Color = 'r';
xl2 = xline(MCmaxFreq/1000);
xl2.Color = 'r';
refl = refline(0,MCminAmp); %minimum amp
refl.Color = 'r';
refl = refline(0,MCmaxAmp); %maximum amp
refl.Color = 'r';

if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.PeakFrequency, Matrixcracks.Amplitude, 'o');
    legend(plus("Hits: ",num2str(HighestIndex)),...
        plus("Matrix cracks: ",...
        num2str(length(Matrixcracks.HitIndex))),...
        'location','south outside');
else
    legend(plus("Hits: ",num2str(HighestIndex)),...
        plus("Matrix cracks: ",...
        num2str(length(Matrixcracks))),'location','south outside');
end
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
    %Debondings
DBSpreadHits = zeros(1,length(time));
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
    %Matrix cracks
MCSpreadHits = zeros(1,length(time));
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
title('Accumulated Acoustic Energy');
%xlim([0 length(TimeVector/Resolution)]);
%ylim([0 max(SpreadEner)+2]);
xlabel('Time [s]');
ylabel('Energy');
plot(time2/Resolution, SpreadEner1);
plot(time2/Resolution, SpreadEner2);
plot(time2/Resolution, SpreadEner3);
legend("Total", "0-200kHz", "200-400kHz", ">400kHz",...
    'location','south outside');
hold off

nexttile %Frequency vs Time vs Amplitude
hold on
PeakFrequencyList = PFreqList/1000;
tbl = table(HitTimeList,PeakFrequencyList,ImpAmpList);
sca = scatter(HitTimeList,PeakFrequencyList,60,ImpAmpList,'.');
sca.CData = tbl.ImpAmpList;
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

nexttile %Duration vs Time
hold on
scatter(HitTimeList, ImpDurList, 60, '.');
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

%Primary 2
figure;
x = linspace(1,100,50);
tiledlayout(2,3);

nexttile %Energy vs Time
hold on
scatter(HitTimeList, ImpEnerList, 60, '.');
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
scatter(ImpPARA1, ImpAmpList, 60, '.');
%refl = refline(0,MCminAmp); %minimum amp
%refl.Color = 'r';
%refl = refline(0,MCmaxAmp); %maximum amp
%refl.Color = 'r';
title('Amplitude vs Load');
%xlim([0 TimeEnd]);
ylim([floor(min(ImpAmpList))*0.9 ceil(max(ImpAmpList))*1.1]);
xlabel('Load [ ]');
ylabel("Amplitude [dB]");
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.Load, Matrixcracks.Amplitude, 'o');
end
hold off

nexttile %Frequency vs Time vs Amplitude
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
title('Frequency vs Time vs Amplitude');
%xlim([0 TimeEnd]);
ylim([0 1000]);
xlabel('Time [s]');
ylabel("Peak Frequency [kHz]");
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.HitTime, Matrixcracks.PeakFrequency, 'o');
end
hold off

nexttile %Energy-Frequency Spectrum
hold on
scatter(PFreqList./1000, ImpEnerList,'.');
xlim([floor(min(PFreqList./1000)/10)*10 ...
    ceil(max(PFreqList./1000)/10)*10]); %nearest 10-number of max sample
ylim([floor(min(ImpEnerList)/10)*10 ...
    ceil(max(ImpEnerList)/10)*10]);
title('Energy-Frequency Spectrum');
xlabel('Peak frequency [kHz]');
ylabel('Energy [dB]');
xl1 = xline(MCminFreq/1000);
xl1.Color = 'r';
xl2 = xline(MCmaxFreq/1000);
xl2.Color = 'r';

if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.PeakFrequency, Matrixcracks.Energy, 'o');
end
hold off

figure
image([0 N], [0 12*10^6], abs(FFTMat));
colorbar

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
    disp("Avg peak frequency:  " + num2str(MCAvgFreq) + "kHz ± " + num2str(MCDevFreq)); 
    disp("Avg peak amplitude:   " + num2str(MCAvgAmp) + "dB ± " + num2str(MCDevAmp)); 
    disp("Avg peak duration:   " + num2str(MCAvgDur) + "μs ± " + num2str(MCDevDur)); 
    disp("Avg peak energy:         " + num2str(MCAvgEner) + " ± " + num2str(MCDevEner)); 
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
    disp("Avg peak frequency:  " + num2str(DBAvgFreq) + "kHz ± " + num2str(DBDevFreq)); 
    disp("Avg peak amplitude:   " + num2str(DBAvgAmp) + "dB ± " + num2str(DBDevAmp)); 
    disp("Avg peak duration:   " + num2str(DBAvgDur) + "μs ± " + num2str(DBDevDur)); 
    disp("Avg peak energy:         " + num2str(DBAvgEner) + " ± " + num2str(DBDevEner)); 
end



