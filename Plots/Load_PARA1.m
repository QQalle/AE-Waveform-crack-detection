figure('name', 'Load and PARA1','Position',[60,60,1400,700])
hold on
% plot(CSVData.Fun_Time, CSVData.Fun_Load,'--')
title('Load and PARA1');
xlim([0 TimeEnd]);
xlabel('Time [s]');
ylabel('Load [N]');
plot(AllValues.HitTime, AllValues.PARA1)
plot(CSVDataOffs.Fun_Time, CSVDataOffs.Fun_Load)
legend('PARA1','offset');
hold off