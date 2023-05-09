figure('name', 'Load vs Time','Position',[60,60,1400,700])
plot(CSVDataOffs.Fun_Time, CSVDataOffs.Fun_Load)
title('Load vs Time');
xlim([0 TimeEnd]);
xlabel('Time [s]');
ylabel('Load [N]');