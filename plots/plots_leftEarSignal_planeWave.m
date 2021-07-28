% SMATBIN: Spherical Microphone Array To Binaural 
%
% Script plots_leftEarSignal_planeWave
% ------------------------------------
%
% Script to generate plots of left ear signals (BRIRs/BRTFs) for the
% minimal working example described in [1] comparing the SMATBIN filter 
% approach with the virtual loudspeaker approach
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
NFFTpw         = 512;           % Length of plane wave applied as test signal impinging on the SMA
soft_limit     = 20;            % Soft-limit for radial filters in dB
sofia_hp       = 1;             % Apply SOFiA's default high-pass at 30 Hz 

% Specify array geometry
radius         = 0.0875;
array_config   = 2;             % 2 for rigid sphear array in SOFiA  
N_grid = [1,3,7];               % Define 3 spatial orders to be calculated and plotted

% Get Lebedev sampling grid for (virtual) SMAs
for k = 1:length(N_grid)
    grid_data_sma{k} = get_sampling_grid('lebedev', N_grid(k));
end

head_orientation = [0 pi/2]; %Frontal head orientation

%% Calculate SMATBIN filters for the different SMA configurations

for k = 1:length(N_grid)
    [smatbin_l{k}, smatbin_r{k}, info{k}] = calc_smatbin_filter(grid_data_sma{k}, N_grid(k), radius, smatbin_filter_length, fs, head_orientation, array_config, soft_limit, sofia_hp);
end

%Internal processing NFFT required for further sample-based comparison of
%SMATBIN and SOFiA
NFFT_proc       = 2*((info{1}.filter_proc_length+NFFTpw));

%% Simulate plane wave for the different SMA configurations

%Incident direction of plane wave
azPW = 0;
colPW = pi/2;

%Circshift parameter of plane wave
ds = (NFFTpw/2)/fs; 

%Plane wave has to be zero-padded to "trick" windowing of decomposed plane waves in sofia_binauralX 
for k = 1:length(N_grid)
    %PW in frequency domain
    [PW{k}, kr_ref{k}] = sofia_swg(radius, grid_data_sma{k}, array_config, fs, NFFTpw, azPW, colPW, 120, ds);
    PW{k} = [PW{k}, conj(PW{k}(:, end-1:-1:2, :))];
    PW{k}(:, NFFTpw/2+1) = real(PW{k}(:, NFFTpw/2+1));
    
    %PW in time domain
    PW_td{k} = real(ifft(PW{k}, [], 2));
    
    %Apply window
    winFkt = hann(size(PW_td{k},2)).';
    for j = 1:size(PW_td{k},1)
        PW_td{k}(j,:) = PW_td{k}(j,:).*winFkt;
    end
    
    %Zero-padded PW in frequency domain
    PW_fd{k} = fft(PW_td{k}, NFFT_proc, 2);
    PW_fd{k} = PW_fd{k}(:, 1:round(end/2)+1);
end

%% Perform binaural rendering with SOFiA (virtual loudspeaker decoding on Lebedev grid)

ncGap = 0; %Default is 16 taps, but 0 required for SMATBIN filter design, so keep similar here.

%Get radial filters for all SMA configurations with NFFT_proc 
fpw  = linspace(0, fs/2, NFFT_proc/2+1);
kpw  = 2*pi*fpw/c;
krpw = kpw*radius;
for k = 1:length(N_grid)
    dn_pw{k} = sofia_mf(N_grid(k), krpw, array_config, soft_limit, 0);
    %Filter kernel truncated to same length as in SMATBIN filter design function
    [dn_pw{k}, dn_pw_kernelSize{k}, dn_pw_Latency{k}] = sofia_rfi(dn_pw{k}, info{1}.filter_proc_length);
end

%Perform binaural rendering for all SMA configurations (Calculate BRIRs)
for k = 1:length(N_grid)
    %SH transform
    PW_nm{k} = sofia_stc(N_grid(k), PW_fd{k}, grid_data_sma{k});
    
    %Binaural decoding using sofia_binauralX (virtual loudspeaker approach)
    [BRTF_sofia_l{k}, BRTF_sofia_r{k}] = sofia_binauralX(PW_nm{k}, dn_pw{k}, head_orientation(1), 1, sofia_hp, 2048, ncGap);
    BRIR_sofia_l{k} = sofia_tdt(BRTF_sofia_l{k});
    BRIR_sofia_r{k} = sofia_tdt(BRTF_sofia_r{k});
end

%% Perform binaural rendering using the SMATBIN filters

NFFTconv = NFFT_proc/2;

%Apply SMATBIN filters for different SMA configurations to SMA signals
for k = 1:length(N_grid)
    smatbin_l_TF{k} = fft(smatbin_l{k}, NFFTconv, 2);
    smatbin_l_TF{k} = smatbin_l_TF{k}(:, 1:round(size(smatbin_l_TF{k},2))/2+1);
    smatbin_r_TF{k} = fft(smatbin_r{k}, NFFTconv, 2);
    smatbin_r_TF{k} = smatbin_r_TF{k}(:, 1:round(size(smatbin_r_TF{k},2))/2+1);
    
    PW_TF{k} = fft(PW_td{k}, NFFTconv, 2);
    PW_TF{k} = PW_TF{k}(:, 1:round(size(PW_TF{k},2))/2+1);
    
    %Convolution in frequency domain + summing
    BRTF_smatbin_l{k} = zeros(1, size(smatbin_l_TF{k},2));
    BRTF_smatbin_r{k} = zeros(1, size(smatbin_r_TF{k},2));
    for chIdx = 1:size(smatbin_l_TF{k},1)
        BRTF_smatbin_l{k} = BRTF_smatbin_l{k} + smatbin_l_TF{k}(chIdx,:) .* PW_TF{k}(chIdx,:);
        BRTF_smatbin_r{k} = BRTF_smatbin_r{k} + smatbin_r_TF{k}(chIdx,:) .* PW_TF{k}(chIdx,:);
    end
    
    %Transform BRTFs to time domain
    BRIR_smatbin_l{k} = sofia_tdt(BRTF_smatbin_l{k});
    BRIR_smatbin_r{k} = sofia_tdt(BRTF_smatbin_r{k});
end

%% For plot: Zero-pad BRIRs to same length

if size(BRIR_sofia_l{1},2) > size(BRIR_smatbin_l{1},2)
    
    for k = 1:length(N_grid)
        zp = zeros(size(BRIR_smatbin_l{k},1), size(BRIR_sofia_l{k},2) - size(BRIR_smatbin_l{k},2));
        BRIR_smatbin_l{k} = [BRIR_smatbin_l{k},zp];
        BRIR_smatbin_r{k} = [BRIR_smatbin_r{k},zp];
    end
    
elseif size(BRIR_sofia_l{1},2) < size(BRIR_smatbin_l{1},2)
    
    for k = 1:length(N_grid)
        zp = zeros(size(BRIR_sofia_l{k},1), size(BRIR_smatbin_l{k},2) - size(BRIR_sofia_l{k},2));
        BRIR_sofia_l{k} = [BRIR_sofia_l{k},zp];
        BRIR_sofia_r{k} = [BRIR_sofia_r{k},zp];
    end
end

%% Cirshift SOFiA BRIRs according to radial_filter_delay
%Copmensation already applied in SMATBIN filter design, i.e., SMATBIN BRIRs do not need to be shifted.

for k = 1:length(N_grid)
    BRIR_sofia_l{k} = circshift(BRIR_sofia_l{k},-info{k}.radial_filter_delay,2);
    BRIR_sofia_r{k} = circshift(BRIR_sofia_r{k},-info{k}.radial_filter_delay,2);
end

%% Save everything

%save data_plots_leftEarSignal

%% Plot parameters

suppMaterial = false;
irPlotWindow = [250 350];
NFFTplot = 4096;

margin1 = [0 0.09];
margin2 = [0.115 0.03];
margin3 = [0.06 0.0075];

gray = [0.7 0.7 0.7];
darkGray = [0.5 0.5 0.5];
blue = [12/255, 94/255, 156/255];
red = [156/255, 27/255, 12/255];
yellow = [249/255, 248/255, 113/255];
cyan = [0/255, 173/255, 162/255];
green = [69/255, 105/255, 31/255];
purple = [75/255, 12/255, 156/255];
lineWidthBack = 0.5;
lineWidthTop = 1.25;
fontSize = 9;

%% Prepare BRIRs for plotting

%Limit amplitude to 1
%Take normalization value from SOFiA as reference so possible amplitude
%variances between SMATBIN and SOFiA are maintained
maxLimAmplitude = 1;
for k = 1:length(N_grid)
    maxAmplitude_sofia(k) = max(abs(BRIR_sofia_l{k}));
    normValue_sofia(k) = maxLimAmplitude/maxAmplitude_sofia(k);
    
    BRIR_sofia_l{k} = BRIR_sofia_l{k} * normValue_sofia(k);
    BRIR_sofia_r{k} = BRIR_sofia_r{k} * normValue_sofia(k);
    
    BRIR_smatbin_l{k} = BRIR_smatbin_l{k} * normValue_sofia(k);
    BRIR_smatbin_r{k} = BRIR_smatbin_r{k} * normValue_sofia(k);
end

%Get spectrum
for k = 1:length(N_grid)
   
    BRIR_sofia_l_Spec{k} = fft(BRIR_sofia_l{k}, NFFTplot, 2);
    BRIR_sofia_l_Spec{k} = BRIR_sofia_l_Spec{k}(:, 1:end/2+1);
    BRIR_sofia_r_Spec{k} = fft(BRIR_sofia_r{k}, NFFTplot, 2);
    BRIR_sofia_r_Spec{k} = BRIR_sofia_r_Spec{k}(:, 1:end/2+1);
    
    BRIR_smatbin_l_Spec{k} = fft(BRIR_smatbin_l{k}, NFFTplot, 2);
    BRIR_smatbin_l_Spec{k} = BRIR_smatbin_l_Spec{k}(:, 1:end/2+1);
    BRIR_smatbin_r_Spec{k} = fft(BRIR_smatbin_r{k}, NFFTplot, 2);
    BRIR_smatbin_r_Spec{k} = BRIR_smatbin_r_Spec{k}(:, 1:end/2+1);
end

%Get frequency vector for plotting
fVec = linspace(0,fs/2,NFFTplot/2+1);

%Get max magnitude values from SOFiA BRIRs as reference
for k = 1:length(N_grid)
   
    maxMagnitude_L(k) = max(20*log10(abs(BRIR_sofia_l_Spec{k})));
    maxMagnitude_R(k) = max(20*log10(abs(BRIR_sofia_r_Spec{k})));
    
    minMagnitude_L(k) = min(20*log10(abs(BRIR_sofia_l_Spec{k})));
    minMagnitude_R(k) = min(20*log10(abs(BRIR_sofia_r_Spec{k})));
   
end
maxMagnitude_L = ceil(max(maxMagnitude_L));
maxMagnitude_R = ceil(max(maxMagnitude_R));
minMagnitude_L = floor(min(minMagnitude_L));
minMagnitude_R = floor(min(minMagnitude_R));

%% Plot for paper

if ~suppMaterial
    close all
    AKf(17,17/2);
    ctn = 0;
    for k = 1:length(N_grid)

        ctn = ctn+1;
        %IRs
        h{ctn} = subtightplot(3,3,ctn,margin1,margin2,margin3);
        p1 = plot(abs(BRIR_sofia_l{k})*-1,'LineStyle','-','LineWidth',lineWidthTop,'Color',red);
        hold on;
        p2 = plot(abs(BRIR_smatbin_l{k}),'LineStyle','-','LineWidth',lineWidthTop,'Color',blue);
        xlim([irPlotWindow(1) irPlotWindow(2)])
        ylim([-maxLimAmplitude,maxLimAmplitude])
        set(gca,'FontSize',fontSize);
        text(irPlotWindow(2)-2,-0.8,['N = ',num2str(N_grid(k))],'FontSize',fontSize,'HorizontalAlignment', 'right');
        if ctn == 1
            set(gca,'yticklabel',{'\pm1','0','1'});
        else
            set(gca,'yticklabel',{'\pm1','0',''});
        end
        if ctn ~= 7
            set(gca,'xticklabel',{[]})
        end
        if ctn == 7
            xlabel('Samples');
            xtickvalues = linspace(irPlotWindow(1), irPlotWindow(2), 3);
            xticks(xtickvalues)
            tickLabels = xtickvalues-irPlotWindow(1);
            set(gca,'xticklabel',{num2str(tickLabels(1)),num2str(tickLabels(2)),num2str(tickLabels(3))});
            set(gca,'yticklabel',{'-1','0',''});
        end

        
        ctn = ctn+1;
        %Magnitude
        h{ctn} = subtightplot(3,3,ctn,margin1,margin2,margin3);
        p1 = semilogx(fVec,20*log10(abs(BRIR_sofia_l_Spec{k})),'LineStyle','-','LineWidth',lineWidthTop,'Color',red);
        hold on;
        p2 = semilogx(fVec,20*log10(abs(BRIR_smatbin_l_Spec{k})),'LineStyle','--','LineWidth',lineWidthTop,'Color',blue);
        p2_legend = semilogx(fVec,20*log10(abs(BRIR_smatbin_l_Spec{k}))+600,'LineStyle','-','LineWidth',lineWidthTop,'Color',blue);
        xlim([20 fs/2])
        ylim([-30,30])
        set(gca,'FontSize',fontSize);
        yticks([-30 -15 0 15 30]);
        if ctn == 2
            set(gca,'yticklabel',{'\pm30','-15','0','15','30'});
            legend([p1,p2_legend],'SOFiA','SMATBIN','Location','NorthWest');
        else
           set(gca,'yticklabel',{'\pm30','-15','0','15',''});
        end
        if ctn ~= 8
            set(gca,'xticklabel',{[]})
        end
        if ctn == 8
            xlabel('Frequency in Hz');
            xticks([10 100 1000 10000])
            set(gca,'xticklabel',{'10','100','1k','10k'});
            set(gca,'yticklabel',{'-30','-15','0','15',''});
        end
        
        ctn = ctn+1;
        %Magnitude Difference
        h{ctn} = subtightplot(3,3,ctn,margin1,margin2,margin3);
        p1 = semilogx(fVec,20*log10(abs(BRIR_sofia_l_Spec{k}))-20*log10(abs(BRIR_smatbin_l_Spec{k})),'LineStyle','-','LineWidth',lineWidthTop,'Color',darkGray);
        xlim([20 fs/2])
        ylim([-1,1])
        set(gca,'FontSize',fontSize);
        yticks([-1 0 1]);
        if ctn == 3
            set(gca,'yticklabel',{'\pm1','0','1'});
        else
           set(gca,'yticklabel',{'\pm1','0',''});
        end
        if ctn ~= 9
            set(gca,'xticklabel',{[]})
        end
        if ctn == 9
            xlabel('Frequency in Hz');
            xticks([10 100 1000 10000])
            set(gca,'xticklabel',{'10','100','1k','10k'});
            set(gca,'yticklabel',{'-1','0',''});
        end
    end

    %yLabel left
    p1 = h{1}.Position;
    p2 = h{7}.Position;
    height=p1(2)+p1(4)-p2(2);
    h3=axes('position',[p2(1)+0.01 p2(2) p2(3) height],'visible','off');
    h_label=ylabel('Amplitude','visible','on','FontSize',fontSize);

    %yLabel middle
    p1 = h{2}.Position;
    p2 = h{8}.Position;
    height=p1(2)+p1(4)-p2(2);
    h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
    h_label=ylabel('Magnitude in dB','visible','on','FontSize',fontSize);
    
    %yLabel right
    p1 = h{3}.Position;
    p2 = h{9}.Position;
    height=p1(2)+p1(4)-p2(2);
    h3=axes('position',[p2(1)+0.01 p2(2) p2(3) height],'visible','off');
    h_label=ylabel('Magnitude difference in dB','visible','on','FontSize',fontSize);

    %Save as PDF
    fileName = ['SMATBINvsSOFiA_LeftEar_HeadAz',num2str(round(rad2deg(head_orientation(1)))),'Col',num2str(round(rad2deg(head_orientation(2)))),'_PWAz',num2str(round(rad2deg(azPW))),'Col',num2str(round(rad2deg(colPW))),'_N_',num2str(N_grid),'_',num2str(smatbin_filter_length),'.pdf'];
    saveas(gcf,fileName);
end

%% Plot for supplementary material

margin1_sm = [0 0.09];
margin2_sm = [0.115 0.09];
margin3_sm = [0.06 0.025];

if suppMaterial
    close all
    AKf(17,17/1.85);
    ctn = 0;
    for k = 1:length(N_grid)

        ctn = ctn+1;
        %IRs
        h{ctn} = subtightplot(3,2,ctn,margin1_sm,margin2_sm,margin3_sm);
        p1 = plot(abs(BRIR_sofia_l{k})*-1,'LineStyle','-','LineWidth',lineWidthTop,'Color',red);
        hold on;
        p2 = plot(abs(BRIR_smatbin_l{k}),'LineStyle','-','LineWidth',lineWidthTop,'Color',blue);
        xlim([irPlotWindow(1) irPlotWindow(2)])
        ylim([-maxLimAmplitude,maxLimAmplitude])
        set(gca,'FontSize',fontSize);
        text(irPlotWindow(2)-2,-0.8,['N = ',num2str(N_grid(k))],'FontSize',fontSize,'HorizontalAlignment', 'right');
        if ctn == 1
            set(gca,'yticklabel',{'\pm1','0','1'});
            text(irPlotWindow(1)+110,1.4,['| HO = (',num2str(round(rad2deg(head_orientation(1)))),'°,',num2str(round(rad2deg(head_orientation(2)))),'°) | K = ',num2str(smatbin_filter_length),' taps |'],'FontWeight','bold','FontSize',fontSize+2,'HorizontalAlignment', 'center');
        else
            set(gca,'yticklabel',{'\pm1','0',''});
        end
        if ctn ~= 5
            set(gca,'xticklabel',{[]})
        end
        if ctn == 5
            xlabel('Samples');
            xtickvalues = linspace(irPlotWindow(1), irPlotWindow(2), 3);
            xticks(xtickvalues)
            tickLabels = xtickvalues-irPlotWindow(1);
            set(gca,'xticklabel',{num2str(tickLabels(1)),num2str(tickLabels(2)),num2str(tickLabels(3))});
            set(gca,'yticklabel',{'-1','0',''});
        end

        
        ctn = ctn+1;
        %Magnitude
        h{ctn} = subtightplot(3,2,ctn,margin1_sm,margin2_sm,margin3_sm);
        p1 = semilogx(fVec,20*log10(abs(BRIR_sofia_l_Spec{k}))-20*log10(abs(BRIR_smatbin_l_Spec{k})),'LineStyle','-','LineWidth',lineWidthTop,'Color',gray);
        hold on;
        p2 = semilogx(fVec,20*log10(abs(BRIR_sofia_l_Spec{k})),'LineStyle','-','LineWidth',lineWidthTop,'Color',red);
        p3 = semilogx(fVec,20*log10(abs(BRIR_smatbin_l_Spec{k})),'LineStyle','--','LineWidth',lineWidthTop,'Color',blue);
        p3_legend = semilogx(fVec,20*log10(abs(BRIR_smatbin_l_Spec{k}))+600,'LineStyle','-','LineWidth',lineWidthTop,'Color',blue);
        xlim([20 fs/2])
        ylim([-30,30])
        set(gca,'FontSize',fontSize);
        yticks([-30 -15 0 15 30]);
        if ctn == 2
            set(gca,'yticklabel',{'\pm30','-15','0','15','30'});
            legend([p2,p3_legend,p1],'SOFiA','SMATBIN','Difference','Location','NorthWest');
        else
           set(gca,'yticklabel',{'\pm30','-15','0','15',''});
        end
        if ctn ~= 6
            set(gca,'xticklabel',{[]})
        end
        if ctn == 6
            xlabel('Frequency in Hz');
            xticks([10 100 1000 10000])
            set(gca,'xticklabel',{'10','100','1k','10k'});
            set(gca,'yticklabel',{'-30','-15','0','15',''});
        end
    end

    %yLabel right
    p1 = h{1}.Position;
    p2 = h{5}.Position;
    height=p1(2)+p1(4)-p2(2);
    h3=axes('position',[p2(1)+0.01 p2(2) p2(3) height],'visible','off');
    h_label=ylabel('Amplitude','visible','on','FontSize',fontSize);

    %yLabel reft
    p1 = h{2}.Position;
    p2 = h{6}.Position;
    height=p1(2)+p1(4)-p2(2);
    h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
    h_label=ylabel('Magnitude in dB','visible','on','FontSize',fontSize);

    %Save as PDF
    fileName = ['SMATBINvsSOFiA_LeftEar_HeadAz',num2str(round(rad2deg(head_orientation(1)))),'Col',num2str(round(rad2deg(head_orientation(2)))),'_PWAz',num2str(round(rad2deg(azPW))),'Col',num2str(round(rad2deg(colPW))),'_N_',num2str(N_grid),'_',num2str(smatbin_filter_length),'_SM.pdf'];
    saveas(gcf,fileName);
end