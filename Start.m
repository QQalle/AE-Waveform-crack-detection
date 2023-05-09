clearvars -except SV && exp
% close all
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
cycle = false; %Use Cycle.m to cycle through multiple experiments
if cycle == false
    experimentNo = '2002'; %Specify which experiment to analyze
else
    experimentNo = SV.experimentNo;
end
run Import_data.m

disable_statistics = true;
disable_read = true;
ApplyHAF = false; %Filter out noise based on FFT integral amplitude

    %Input hardware calibrations
PT = 20*10^-6; %Pre-trigger
PDT = 35; %Peak Definition Time
HDT = 150; %Hit Definition Time
HLT = 300; %Hit Lockout Time
Fs = 5*10^6; %Sample frequency (Hz)

    %Input software parameters
Total = length(ASCIIOutPut.data)/Fs;
TimeEnd = 90; %Experiment cutoff time [s]
SampleNumber = 1; %Matrix crack no. to sample !OBS ERROR IF > TOTAL!
HAFfilter = -1500; %Recommended: -3000 - 100 Default: -1500

    %Calibrate matrix crack definition
MCminFreq = 400*10^3; %[Hz]
MCmaxFreq = 1000*10^3; %[Hz]
MCminAmp = 0; %[dB]
MCmaxAmp = 200; %[dB]
MCminEner = 0; %[kJ]
MCmaxEner = 10^16; %[aJ]
MCminDur = 0; %[μs]
MCmaxDur = 2000;
MCminCount = 0;
MCmaxCount = 50;
MCminRise = 0; %[s]
MCmaxRise = 45; %[s]
MCstr = 3; % [%] What is expected start strain for matrix cracks
EnergyCap = 1*10^7; %To exclude anomalies

    %Calibrate debonding definition
DBminFreq = 240*10^3; %[Hz]
DBmaxFreq = 310*10^3; %[Hz]
DBminAmp = 45; %[dB]
DBmaxAmp = 65; %[dB]
DBminEner = 0; %[kJ]
DBmaxEner = 10^12; %[aJ]
DBminDur = 0; %[μs]
DBminCount = 100;
DBminRise = 0; %[s]

CheckVariable = "duration";
%Varibles: "peak frequency", "amplitude", "duration", "energy",
%          "counts," "rise time", "parametric 1"
CheckRangeMIN = 4000;
CheckRangeMAX = 400000;


run Calculations.m %Load waveforms and calculations

if cycle == false %Make graphs
    %run Plots.m 
    %%% Other %%%
%     run Plots/Amplitude_PARA1.m
%     run Plots/Sample_Waveform.m
    %%% Time graphs %%%
%     run Plots/Duration_Time.m
%     run Plots/Energy_Time.m
%     run Plots/Load_Time.m
%     run Plots/Stress_Time.m
%     run Plots/Strain_Time.m
%     run Plots/Counts_Time.m
%     run Plots/Risetime_Time.m
    %%% Frequency %%%
%     run Plots/Amplitude_Frequency.m
%     run Plots/Frequency_Time_Amplitude.m
%     run Plots/Frequency_Duration_Energy.m
%     run Plots/Energy_Frequency.m
    %%% Advanced %%%
    run Plots/Spectrogram.m %(HAF)
    run Plots/Hitcounter.m
    run Plots/Cumulative_Energy.m
    run Plots/Cumulative_Energy_Stress.m
    run Plots/Stacked_Hits.m
    run Plots/Energy_Derivative.m
    run Plots/Stress_Hits.m
    %%% Debug %%%
%     run Plots/Load_PARA1.m


end

if disable_statistics == false
    run Statistics.m %Compute statistics of hits
end
disp('it worked! =)')