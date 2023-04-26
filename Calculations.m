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

PARAStart = 1;
while ImpPARA1(PARAStart) <= 100/10000 %Check for index when test starts
    PARAStart = PARAStart + 1;
end

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
        FunTime = TimeVector;
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
time2 = (0:length(TimeVector)-1);
TestStart = HitTimeList(PARAStart);

CSVData.Fun_Time = CSVData.Fun_Time + TestStart;
CSVDataOffs = array2table(NaN(length(TimeVector),width(CSVData)),...
    'VariableNames',CSVData.Properties.VariableNames);
CSVDataOffs.Fun_Time = transpose(0:length(time2)-1)/Resolution;
for i = 1 : length(CSVData.Fun_Time)
    t = round(CSVData.Fun_Time(i),1)+1;
    index = t*Resolution;
    CSVDataOffs(index,2:end) = CSVData(i,2:end);
       
end
CSVDataOffs = fillmissing(CSVDataOffs,"next");
%CSVDataOffs = fillmissing(CSVDataOffs,"previous");
PullStop = CSVData.Fun_Time(end) + TestStart; %(s) When tensile test ends (stops pulling)
disp("Tensile test end =" + num2str(PullStop) + "s");

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