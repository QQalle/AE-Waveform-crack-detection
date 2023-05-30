figure('name', 'Energy-time derivative','Position',[60,60,1400,700])
disp("calculating Energy-time derivative...")
hold on
plot(time2/Resolution, DerivEner);
% plot(time2/Resolution, DerivEner1,'--');
% plot(time2/Resolution, DerivEner2,'--');
% plot(time2/Resolution, DerivEner3,'--');
title('Energy Derivative vs Time');
xlabel('Time [s]');
ylabel('Energy derivative');
% if istable(Matrixcracks) %If true = there are matrixcracks
%     for i = 1 : length(Matrixcracks.HitTime)
%         xli = xline(Matrixcracks.HitTime(i));
%         xli.Color = [1 0 0];
%     end
% end
% findpeaks(DerivEner,Fs,'MinPeakHeight',500)
disp(append('Peaks: ',num2str(EnerPeaksNo)));
disp(append('Peaks divided: ',num2str(EnerPeaksNo2)));
hold off