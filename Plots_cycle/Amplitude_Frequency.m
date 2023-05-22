figure('name', 'Amplitude-Frequency','Position',[60,1,1000,1000])
% t.TileSpacing = 'none';
hold on
for k = 1 : exp
%     nexttile
    figure('name', 'Amplitude-Frequency','Position',[60,1,1000,1000])
    plot(SV.PeakFrequency{k}./1000, SV.Amplitude{k}, '.');
%     r = refline(0,45);
%     r.Color = 'r';
    title(append(num2str(Order2(k)),': Amp-Freq, ',num2str(round(SV.Table{k,1})),' MPa'));
    xlim([0, 1000]);
    xticks(0:200:1000)
%     xlabel('Peak frequency [kHz]');
%     ylabel('Amplitude [dB]');
end
