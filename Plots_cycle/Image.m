figure('name', 'Frequency vs Stacked hits','Position',[60,60,1400,700])
disp("Stacking hits...");
hold on
multiplier = 1;
for k = 1 : exp
%     nexttile
    Spreadf = cell(1,multiplier*ceil(max(SV.PeakFrequency{k}./1000)));
    Spreadf2 = [];
    PeakFrequency = flip(SV.PeakFrequency{k});
    for i = 1 : length(PeakFrequency)
        Spreadf{1,round(PeakFrequency(i)./1000*multiplier)} = ...
            [Spreadf{1,round(PeakFrequency(i)./1000*multiplier)};...
            i*SV.TimeEnd(k)/length(SV.Experiments)];
    end
    % Spreadf = Spreadf(1,1:(1000*multiplier)); %Cap the frequency to 1MHz
    for i = 1 : length(Spreadf)
        len = length(Spreadf{1,i});
        if len ~= 0
            if len == 1
                Spreadf2(1,i) = 1;
            else
                Spreadf2(1:len,i) = Spreadf{1,i};
            end
        end
    end
    axis([0 1000 0.5 inf]);
    img = imagesc(Spreadf2/multiplier);
    cc = colorbar;
    title(cc,'Time [s]')
    shading interp
    mycolormap = colormap(hsv);
    % mycolormap(1, :) = [1,1,1]; %White background
    mycolormap(1, :) = [0, 0, 0]; %Black background
    % bound = MCstrboundtime*height(mycolormap)/TimeEnd;
    % mycolormap(1:round(bound), :) = zeros(round(bound), 3); %Reduced lower values
    colormap(mycolormap);
    title('Frequency vs Stacked hits');
    xlabel('Peak frequency [kHz]');
    ylabel('Hits');
    % set(gca,'Color','k')
%     set(gca,'FontSize',14)
end
hold off
