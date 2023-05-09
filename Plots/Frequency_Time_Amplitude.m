figure('name', 'Frequency vs Time vs Amplitude','Position',[60,60,1400,700])
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