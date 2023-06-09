fclose all
clear
warning('off','all')
warning

experimentNo = [];
% Experiments = ["2001","2002","2003","2004","2005",...
%    "2006","2007","2008","2009","3001","3002","3003","3004","3005"];
% Experiments = ["2005","2006","2004","2007","2003","2002","3003",...
%      "2001","3001"]; %Sorted one of each stress level
% Experiments = ["2005","2006","2004","2007","2003","2002","3003",...
%      "3004","2001","3002","3001","3005","2009","2008"]; %falling order
% Experiments = ["3003","3004","3002","3001","3005","2008","2009"]; %Crack amount order
                
% Experiments = ["3001", "3002","3003","3004"];

% For the Crack_energy graphs
Experiments = ["3003", "3004", "3002", "3001", "3005", "2009", "2008"]; 

SV = struct('experimentNo',experimentNo,'Experiments',Experiments);
for exp = 1 : length(SV.Experiments)
        %Variables to update
    experimentNo = SV.Experiments(exp);
    SV.experimentNo = experimentNo; %Variables to save
%     disp(experimentNo)
    if exp == 1
        StressThres = [];
        TEnerThres = [];
        SpreadEner = [];
        MPa = [];
        Table = table(MPa,StressThres,TEnerThres,SpreadEner);
        SV.Table = Table;
        SV.time2 = cell(1,length(SV.Experiments));
        SV.SpreadEner = cell(1,length(SV.Experiments));
        SV.SpreadEner1 = cell(1,length(SV.Experiments));
        SV.SpreadEner2 = cell(1,length(SV.Experiments));
        SV.SpreadEner3 = cell(1,length(SV.Experiments));
        SV.EnergyTable = table();
        SV.EnergyTable1 = table();
        SV.EnergyTable2 = table();
        SV.EnergyTable3 = table();
        SV.TimeEnd = [];
        SV.Fun_Time = cell(1,length(SV.Experiments));
        SV.Fun_TensileStress = cell(1,length(SV.Experiments));
        SV.Resolution = [];
        SV.HitTimeList = [];
        SV.SpreadEnerAPE = [];
        SV.SpreadEnerAPETime = [];
        SV.PullStop = [];
        SV.PullStopIndex = [];
        SV.PullStopIndexPlus = [];
        SV.PullStopIndexStress = [];
        SV.EstCracks = table();
        SV.EstCracksTotal = cell(1,length(SV.Experiments));
        SV.PeakFrequency = cell(1,length(SV.Experiments));
        SV.Amplitude = cell(1,length(SV.Experiments));
        SV.DerivEner = cell(1,length(SV.Experiments));
    end
    run Start.m
    SV.Table.MPa(exp) = max(CSVDataOffs.Fun_TensileStress);
    if isempty(StressThres) == 0
        SV.Table.StressThres(exp) = StressThres;
        SV.Table.TEnerThres(exp) = TEnerThres;
        SV.Table.SpreadEner(exp) = max(SpreadEner);
    end
    SV.time2{exp} = time2;
    SV.SpreadEner{exp} = SpreadEner;
    SV.SpreadEner1{exp} = SpreadEner1;
    SV.SpreadEner2{exp} = SpreadEner2;
    SV.SpreadEner3{exp} = SpreadEner3;
    SV.TimeEnd(exp) = TimeEnd;
    SV.Fun_Time{exp} = CSVDataOffs.Fun_Time;
    SV.Fun_TensileStress{exp} = CSVDataOffs.Fun_TensileStress;
    SV.Resolution(exp) = Resolution;
    SV.HitTimeList{exp} = HitTimeList;
    SV.SpreadEnerAPE{exp} = SpreadEnerAPE;
    SV.SpreadEnerAPETime{exp} = SpreadEnerAPETime;
    SV.PullStop(exp) = PullStop;
    SV.PullStopIndex(exp) = find(SV.time2{exp}/Resolution >= PullStop,1);
    if length(SpreadEner) >= SV.PullStopIndex(exp) + Resolution * 15
        SV.PullStopIndexPlus(exp) = SV.PullStopIndex(exp) + Resolution * 10; %15 for 15 sec
    else
        SV.PullStopIndexPlus(exp) = length(SpreadEner);
    end
    SV.PullStopIndexStress(exp) = find(SV.Fun_Time{exp} >= SV.PullStop(exp), 1);
    SV.EnergyTable{exp,1} = experimentNo;
    SV.EnergyTable{exp,2} = max(SV.SpreadEner{exp});
    SV.EnergyTable{exp,3} = SV.SpreadEner{exp}(SV.PullStopIndex(exp));
    SV.EnergyTable{exp,4} = SV.SpreadEner{exp}(SV.PullStopIndexPlus(exp));
    SV.EnergyTable1{exp,1} = experimentNo;
    SV.EnergyTable1{exp,2} = max(SV.SpreadEner1{exp});
    SV.EnergyTable1{exp,3} = SV.SpreadEner1{exp}(SV.PullStopIndex(exp));
    SV.EnergyTable1{exp,4} = SV.SpreadEner1{exp}(SV.PullStopIndexPlus(exp));
    SV.EnergyTable2{exp,1} = experimentNo;
    SV.EnergyTable2{exp,2} = max(SV.SpreadEner2{exp});
    SV.EnergyTable2{exp,3} = SV.SpreadEner2{exp}(SV.PullStopIndex(exp));
    SV.EnergyTable2{exp,4} = SV.SpreadEner2{exp}(SV.PullStopIndexPlus(exp));
    SV.EnergyTable3{exp,1} = experimentNo;
    SV.EnergyTable3{exp,2} = max(SV.SpreadEner3{exp});
    SV.EnergyTable3{exp,3} = SV.SpreadEner3{exp}(SV.PullStopIndex(exp));
    SV.EnergyTable3{exp,4} = SV.SpreadEner{exp}(SV.PullStopIndexPlus(exp));
    SV.EstCracks{exp,1} = experimentNo;
    SV.EstCracks{exp,2} = round((2*10^-6)*SV.EnergyTable{exp,3}-10.114);
    disp(append('Estimated cracks: ',num2str(SV.EstCracks{exp,2})));
    SV.EstCracks{exp,3} = 0.0083*SV.Table.MPa(exp)^2-3.7326...
        *SV.Table.MPa(exp)+414.22;
    SV.EstCracks{exp,4} = mean([SV.EstCracks{exp,2},SV.EstCracks{exp,3}]);
    for i = 1:length(SV.SpreadEner{exp})
        if round((2*10^-6)*SV.SpreadEner{exp}(i)-10.114) > 0
            SV.EstCracksTotal{exp}(i) = ...
                round((2*10^-6)*SV.SpreadEner{exp}(i)-10.114);
        else
            SV.EstCracksTotal{exp}(i) = 0;
        end
    end
    SV.PeakFrequency{exp} = PFreqList;
    SV.Amplitude{exp} = HAFImpAmpList;
    SV.DerivEner{exp} = DerivEner;
end
%%
Markers = ["o-","*-","x-","square-","diamond-","^-","v-","<-",...
    ">-","pentagram-","hexagram-"];
Markers = [Markers Markers Markers Markers];
% Order = varfun(@(x) (num2str(x)),flip((1:1:length(SV.Experiments))));
Order = flip((1:1:length(SV.Experiments)));
Order2 = [14,13,12,11,10,9,8,6,4];
Legendtext = plus(num2str(Order'),...
    plus(": ",...
    plus(num2str(round(SV.Table.MPa)),...
    "MPa")));
% Legendtext2 = plus(num2str(Order2'),...
%     plus(": ",...
%     plus(num2str(round(SV.Table.MPa)),...
%     "MPa")));
capped = true;
if capped == true
    captext = ', filtered';
else
    captext = '';
end

SV.EnergyTable{:,1} = Legendtext;
SV.EnergyTable1{:,1} = Legendtext;
SV.EnergyTable2{:,1} = Legendtext;
SV.EnergyTable3{:,1} = Legendtext;
SV.EstCracks{:,1} = Legendtext;

figure('name', 'Cumulative Acoustic Energy vs Stress','Position',...
    [60,60,1400,700])
hold on
for k = 1 : exp
%     yyaxis left
%     plot(SV.time2{k}/SV.Resolution(k), SV.SpreadEner{k});
%     PullStopInd = find(SV.Fun_Time{k} >= SV.PullStop(k), 1);
    PullStopInd = find(SV.Fun_Time{k} >= SV.PullStop(k), 1);
    p = plot(SV.Fun_TensileStress{k}, SV.SpreadEner{k},Markers(k),...
        'MarkerIndices',PullStopInd-2);
%     yyaxis right
%     plot(SV.Fun_Time{k}, SV.Fun_TensileStress{k})
end
title(append('Cumulative Energy vs Stress',captext));
% xlim([0 max(SV.Fun_Time)]);
% xlabel('Time [s]');
xlabel('Stress [MPa]');
ylabel('Energy [aJ]');
legend(Legendtext,'location','east outside');
grid on
yticks((0:1:15)*10^7)
xticks(0:50:800)
set(gca,'FontSize',14)
hold off

% Accumulated Energy after PullStop v. Time
figure('name', 'Proportion of total accumulated Energi after pullstop vs Time',...
    'Position',[60,60,1400,700])
hold on
for k = 1 : exp
    plot(SV.SpreadEnerAPETime{k}, 100*SV.SpreadEnerAPE{k},Markers(k),...
        'MarkerIndices',length(SV.SpreadEnerAPETime{k}));
end
title(append('Proportion of Total Accumulated Energy After Pullstop vs Time',...
    captext));
% xlim([0 max(SV.Fun_Time)]);
% xlabel('Time [s]');
xlabel('Time [sec]');
ylabel('Percentage of total Energy after Pullstop');
legend(Legendtext,'location','east outside');
grid on
ytickformat('percentage')
set(gca,'FontSize',14)
hold off

% Stress v. Hits
figure('name','Hits vs Stress','Position',[60,60,1400,700])
hold on
for k = 1 : exp
    SpreadHits2 = zeros(1,length(SV.time2{k}));
    for i = 1 : length(SV.HitTimeList{k})
        SpreadHits2(ceil(SV.HitTimeList{k}(i))*SV.Resolution(k):end) = ...
            SpreadHits2(ceil(SV.HitTimeList{k}(i))*SV.Resolution(k):end)+1;
    end
    SmoothSpreadHits2 = smooth(SpreadHits2,100);
    PullStopInd = find(SV.Fun_Time{k} >= SV.PullStop(k), 1);
    plot(SV.Fun_TensileStress{k}, SmoothSpreadHits2,Markers(k),...
        'MarkerIndices',PullStopInd-2);
end
title(append('Hits vs Stress',captext));
xlabel('Stress [MPa]');
ylabel('Hits');
legend(Legendtext,'location','east outside');
grid on
xticks(0:50:800)
set(gca,'FontSize',14)
hold off

% Hits v. time (after pullstop)
figure('name','Proportion of Total Hits After Pullstop vs Time','Position',[60,60,1400,700])
hold on
for k = 1 : exp
    SpreadHits2 = zeros(1,length(SV.time2{k}));
    for i = 1 : length(SV.HitTimeList{k})
        SpreadHits2(ceil(SV.HitTimeList{k}(i))*SV.Resolution(k):end) = ...
            SpreadHits2(ceil(SV.HitTimeList{k}(i))*SV.Resolution(k):end)+1;
    end
    SmoothSpreadHits2 = smooth(SpreadHits2,100);
    PullStopInd = find(SV.Fun_Time{k} >= SV.PullStop(k), 1);
    SmoothSpreadHits2origo = (SmoothSpreadHits2(PullStopInd:end)...
        - min(SmoothSpreadHits2(PullStopInd:end)))...
        ./ max(SmoothSpreadHits2);
    plot(SV.Fun_Time{k}(1:length(SmoothSpreadHits2origo))...
        , 100*SmoothSpreadHits2origo,Markers(k),...
        'MarkerIndices',length(SV.SpreadEnerAPETime{k}));
end
title(append('Proportion of Total Hits After Pullstop vs Time',captext));
xlabel('Time [S]');
ylabel('Percentage of Hits after Pullstop');
ytickformat('percentage')
legend(Legendtext,'location','east outside');
grid on
% yticks((0:1:15)*10^7)
% xticks(0:50:800)
set(gca,'FontSize',14)
hold off

% run Plots_cycle/A_E_Stress_0_200kHz.m
% run Plots_cycle/A_E_Stress_200_400kHz.m
% run Plots_cycle/A_E_Stress_400_infkHz.m
% run Plots_cycle/Estimate_cracks.m
% run Plots_cycle/Amplitude_Frequency.m
run Plots_cycle/DerivEnerMEAN.m

disp('it worked! =)')