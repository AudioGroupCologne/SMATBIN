% SMATBIN: Spherical Microphone Array To Binaural 
%
% Script plots_SMATBINfilter
% --------------------------
%
% Script to generate plots of left ear SMATBIN filters as presented in [1]
% 
% Dependencies:
% -------------
% SOFiA toolbox: https://github.com/AudioGroupCologne/SOFiA
% AKtools: https://www.ak.tu-berlin.de/menue/publications/open_research_tools/aktools/
%
% References:
% -------------
% [1] J.M. Arend, T. Lübeck, and C. Pörschmann, 
% "Efficient binaural rendering of spherical microphone array data by linear filtering", 
% Submitted for publication 
%
% (C) 2021 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing 
%%
clear all; close all; clc

% Define array configuartions:
c              = 343;    
fs             = 48000;     
smatbin_filter_length = 2048;   % Length of SMATBIN filters    
soft_limit     = 20;            % Soft-limit for radial filters in dB
sofia_hp       = 1;             % Apply SOFiA's default high-pass at 30 Hz 

% Specify array geometry
radius         = 0.0875;
array_config   = 2;             % 2 for rigid sphear array in SOFiA  
N_grid         = 1;             % Spatial order of SMA

% Get Lebedev sampling grid for (virtual) SMA
grid_data_sma = get_sampling_grid('lebedev', N_grid);

%Define head orientations (here only 4 horizontal head orientations)
head_orientations = get_sampling_grid('horizontal');
head_orientations = head_orientations(1:90:end, :);

%% Calculate SMATBIN filters for the SMA configuration

[smatbin_l, smatbin_r, info] = calc_smatbin_filter(grid_data_sma, N_grid, radius, smatbin_filter_length, fs, head_orientations, array_config, soft_limit, sofia_hp);

%% Plot parameters

gridDeg = grid_data_sma*180/pi; %Transform grid to degree
irCut = 100; %Define window size for plotting (IR is not getting truncated)
NFFTplot = 4096;
if NFFTplot < smatbin_filter_length
    NFFTplot = smatbin_filter_length;
end
nPoints = length(grid_data_sma);
headID = 1; %Define head orientation to be plotted. Here 1 = frontal head orientation

margin1 = [0 0.09];
margin2 = [0.05 0.01];
margin3 = [0.06 0.025];

gray = [0.7 0.7 0.7];
blue = [12/255, 94/255, 156/255];
red = [156/255, 27/255, 12/255];
lineWidthBack = 0.5;
lineWidthTop = 1.25;
fontSize = 8;
fontSizeBoxText = 8;

%% Prepare SMATBIN filter for plotting

%Limit amplitude to 1
maxLimAmplitude = 1;
maxAmplitude = max(max(max(smatbin_l(:,:,:)))); %Is the same for any head id..
smatbin_l = smatbin_l * maxLimAmplitude/maxAmplitude;

%Get spectrum
bin_sma_fir_TF_l = fft(smatbin_l, NFFTplot, 2);
bin_sma_fir_TF_l = bin_sma_fir_TF_l(:, 1:round(size(bin_sma_fir_TF_l, 2)/2 +1), :);

%Get respective frequency vector
fVec = linspace(0,fs/2,NFFTplot/2+1);

%Get max and min magnitude
maxMagnitude = round(max(max(max(20*log10(abs(bin_sma_fir_TF_l))))));
minMagnitude = floor(min(min(min(20*log10(abs(bin_sma_fir_TF_l))))));

%Calculate groupdelay for all filters
clear phiAll gdAll;
for n = 1:nPoints
        phiAll(n,:) = unwrap(angle(bin_sma_fir_TF_l(n,:,headID)));
        gdAll(n,:) = -diff(phiAll(n,:).') / (2*pi) * NFFTplot/(fs/2);
end

%Get respective frequency vector
fVecGD = linspace(0,fs/2,NFFTplot/2);

%% Plot for paper

close all
fig = AKf(17,17);
ctn = 0;
for n = 1:nPoints

    ctn = ctn+1;
    %IRs
    h{ctn} = subtightplot(6,3,ctn,margin1,margin2,margin3);
    plot(squeeze(smatbin_l(:,1:irCut,headID).'),'Color',gray,'LineWidth',lineWidthBack);
    hold on;
    plot(squeeze(smatbin_l(n,1:irCut,headID)),'LineWidth',lineWidthTop,'Color',blue);
    xlim([0 irCut])
    ylim([-maxLimAmplitude,maxLimAmplitude])
    set(gca,'FontSize',fontSize);
    text(irCut-3,0.8,['Az = ',num2str(gridDeg(n,1)),'°, Col = ',num2str(gridDeg(n,2)),'°'],'FontSize',fontSize,'HorizontalAlignment', 'right');
    if ctn == 1
        set(gca,'yticklabel',{'\pm1','-0.5','0','0.5','1'});
    else
        set(gca,'yticklabel',{'\pm1','-0.5','0','0.5',''});
    end
    if ctn ~= 16
        set(gca,'xticklabel',{[]})
    end
    if ctn == 16
        xlabel('Samples');
        set(gca,'yticklabel',{'-1','-0.5','0','0.5',''});
    end


    ctn = ctn+1;
    %Magnitude
    h{ctn} = subtightplot(6,3,ctn,margin1,margin2,margin3);
    semilogx(fVec,20*log10(abs(bin_sma_fir_TF_l(:,:,headID))),'Color',gray,'LineWidth',lineWidthBack);
    hold on;
    semilogx(fVec,20*log10(abs(bin_sma_fir_TF_l(n,:,headID))),'LineWidth',lineWidthTop,'Color',blue);
    xlim([0 fs/2])
    ylim([-30,30])
    set(gca,'FontSize',fontSize);
    yticks([-30 -15 0 15 30]);
    if ctn == 2
        set(gca,'yticklabel',{'\pm30','-15','0','15','30'});
    else
       set(gca,'yticklabel',{'\pm30','-15','0','15',''});
    end
    if ctn ~= 17
        set(gca,'xticklabel',{[]})
    end
    if ctn == 17
        xlabel('Frequency in Hz');
        xticks([100 1000 10000])
        set(gca,'xticklabel',{'100','1k','10k'});
        set(gca,'yticklabel',{'-30','-15','0','15',''});
    end

    
    ctn = ctn+1;
    %Group Delay
    h{ctn} = subtightplot(6,3,ctn,margin1,margin2,margin3);
    semilogx(fVecGD,gdAll*1000,'Color',gray,'LineWidth',lineWidthBack);
    hold on;
    semilogx(fVecGD,gdAll(n,:)*1000,'LineWidth',lineWidthTop,'Color',blue);
    xlim([0 fs/2])
    ylim([-15,15])
    set(gca,'FontSize',fontSize);
    yticks([-15 -10 -5 0 5 10 15]);
    if ctn == 3
        set(gca,'yticklabel',{'\pm15','-10','-5','0','5','10','15'});
    else
       set(gca,'yticklabel',{'\pm15','-10','-5','0','5','10',''});
    end
    if ctn ~= 18
        set(gca,'xticklabel',{[]})
    end
    if ctn == 18
        xlabel('Frequency in Hz');
        xticks([100 1000 10000])
        set(gca,'xticklabel',{'100','1k','10k'});
        set(gca,'yticklabel',{'-15','-10','-5','0','5','10',''});
    end
end

%yLabel left
p1 = h{1}.Position;
p2 = h{16}.Position;
height=p1(2)+p1(4)-p2(2);
h3=axes('position',[p2(1)+0.01 p2(2) p2(3) height],'visible','off');
h_label=ylabel('Amplitude','visible','on','FontSize',fontSize);

%yLabel middle
p1 = h{2}.Position;
p2 = h{17}.Position;
height=p1(2)+p1(4)-p2(2);
h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
h_label=ylabel('Magnitude in dB','visible','on','FontSize',fontSize);

%yLabel right
p1 = h{3}.Position;
p2 = h{18}.Position;
height=p1(2)+p1(4)-p2(2);
h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
h_label=ylabel('Group delay in ms','visible','on','FontSize',fontSize);

%Save as PDF
fileName = ['SMATBIN_LebedevN1_HeadAz0Col90_GD_',num2str(smatbin_filter_length),'.pdf'];
saveas(gcf,fileName);
