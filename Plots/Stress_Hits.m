figure('name','Stress-Hits','Position',[60,60,1400,700])
SpreadHits2 = zeros(1,length(time2));
for i = 1 : length(HitTimeList)
    SpreadHits2(ceil(HitTimeList(i))*Resolution:end) = ...
        SpreadHits2(ceil(HitTimeList(i))*Resolution:end)+1;
end
SmoothSpreadHits2 = smooth(SpreadHits2,100);
plot(CSVDataOffs.Fun_TensileStress, SmoothSpreadHits2);
title('Stress-Hits');
xlabel('Stress [MPa]');
ylabel('Hits');
% plot(time2/Resolution, SmoothSpreadHits2);
% plot(time2/Resolution, CSVDataOffs.Fun_Load);