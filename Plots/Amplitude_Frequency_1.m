figure('name', 'Frequency_Duration_Energy','Position',[60,60,1400,700])
hold on
scatter(PFreqList./1000, HAFImpDurList, 30, HAFImpEnerList, '.');
colorbar
% xlim([floor(min(PFreqList./1000)/10)*10 ...
%     ceil(max(PFreqList./1000)/10)*10]); %nearest 10-number of max sample
% ylim([floor(min(HAFImpAmpList)/10)*10 ...
%     ceil(max(HAFImpAmpList)/10)*10]);
title('Frequency_Duration_Energy');
xlabel('Peak frequency [kHz]');
ylabel('Duration [μs]');