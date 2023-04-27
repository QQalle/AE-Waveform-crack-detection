clear
close all
%{
    Instructions
Choose no less than 2000 data collection length
Make sure to format files as following example:
AEwin Output export: "EXP[experiment number].txt"
AEwin Waveforms folder: "EXP[experiment number]"
Tensile machine output CSV: "Specimen_RawData_[experiment number].csv"

    1001:
Fs = 5*10^6 %(Hz)
PullStop = 56.382 %(s)
    1002:
Fs = 5*10^6 %(Hz)
PullStop = 51.838 %(s)
    1003:
Fs = 10*10^6 %(Hz)
PullStop = 48.27 %(s)
TimeEnd = 124
%}

experimentNo = '1003'; %Specify which experiment to analyze
run Import_data.m

    %Input hardware calibrations
PT = 20*10^-6; %Pre-trigger
PDT = 35; %Peak Definition Time
HDT = 150; %Hit Definition Time
HLT = 300; %Hit Lockout Time
Fs = 10*10^6; %Sample frequency (Hz)

    %Input software parameters
Total = length(ASCIIOutPut.data)/Fs;
TimeEnd = 80; %Experiment cutoff time [s]
SampleNumber = 1; %Matrix crack no. to sample !OBS ERROR IF > TOTAL!
ApplyHAF = false; %Filter out noise based on FFT integral amplitude
HAFfilter = -1500; %Recommended: -3000 - 100 Default: -1500

    %Calibrate matrix crack definition
MCminFreq = 75*10^3; %[Hz]
MCmaxFreq = 180*10^3; %[Hz]
MCminAmp = 60; %[dB]
MCmaxAmp = 99; %[dB]
MCminEner = 0; %[kJ]
MCmaxEner = 10^12; %[aJ]
MCminDur = 1500; %[μs]

    %Calibrate debonding definition
DBminFreq = 240*10^3; %[Hz]
DBmaxFreq = 310*10^3; %[Hz]
DBminAmp = 45; %[dB]
DBmaxAmp = 65; %[dB]
DBminEner = 0; %[kJ]
DBmaxEner = 10^12; %[aJ]
DBminDur = 0; %[μs]

run Calculations.m
run Plots.m