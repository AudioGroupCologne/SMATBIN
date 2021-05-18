% SMATBIN: Spherical Microphone Array To Binaural 
%
% Script DEMO1_plane_wave_sma
% ---------------------------
%
% DEMO 1: Comparison of binaural rendering of a single plane wave impinging 
% on a simulated spherical microphone array using SMATBIN filters and the 
% SOFiA rendering chain 
%
% References:
% -------------
% [1] J.M. Arend, T. Lübeck, and C. Pörschmann, 
% "Efficient binaural rendering of spherical microphone array data by linear filtering", 
% Submitted for publication 
%
% (C) 2021 by TL, Tim Lübeck
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing 
%%
clear all; close all; clc
export_plots = 1;

% Define array configuartions:
c              = 343;    
fs             = 48000;     
smatbin_filter_length = 2048;   % Length of SMATBIN filters    
NFFT_pw         = 512;          % Length of plane wave applied as test signal impinging on the SMA
soft_limit     = 40;            % Soft-limit for radial filters in dB
sofia_hp       = 1;             % Apply SOFiA's default high-pass at 30 Hz 

% Specify array geometry
radius          = 0.0875;
array_body      = 2;            % 2 for rigid sphear array in SOFiA      
N_grid          = 1;            % Spatial order of SMA

% Get Lebedev sampling grid for (virtual) SMA
[grid_data_sma, ~, N_grid, ~] = get_sampling_grid('lebedev', N_grid);

% Define head orientations to be rendered
head_orientations = get_sampling_grid('horizontal');
head_orientations = head_orientations(1:90:end, :);
%head_orientations = [0, 0];

%% Calculate SMATBIN filters for the SMA configuration

[smatbin_l, smatbin_r, info] = calc_smatbin_filter(grid_data_sma, N_grid, radius, ...
                                        smatbin_filter_length, fs, head_orientations, array_body, soft_limit, sofia_hp);

% calculate processing NFFT for conventiononal binaural rendering 
NFFT_proc = 2*((info.filter_proc_length + NFFT_pw));

%% Simulate plane wave incident on the SMA

% Incident direction of plane wave
az   = 0*(pi/180); 
col  = 90*(pi/180);   

% Circshift parameter of plane wave
ds = (NFFT_pw/2)/fs;      

[DRTFs_simulated, ~] = sofia_swg(radius, grid_data_sma, array_body, fs, ...
                             NFFT_pw, az, col, 120, ds);
                         
DRTFs_simulated = [DRTFs_simulated, conj(DRTFs_simulated(:, size(DRTFs_simulated, 2)-1:-1:2, :))];
DRTFs_simulated(:, NFFT_pw/2+1) = real(DRTFs_simulated(:, NFFT_pw/2+1));
drirs_simulated = real(ifft(DRTFs_simulated, [], 2));

% window
winFkt = hann(size(drirs_simulated, 2)).';
for j = 1:size(drirs_simulated, 1)
    drirs_simulated(j, :) = drirs_simulated(j, :) .* winFkt;
end
    
%% Evaluation parameters

head_or_eval = [0, pi/2]; % Head orientation for evaluation
head_or_idx = get_nn_idx(head_orientations, head_or_eval(1), head_or_eval(2)); % idx of head orientation for evaluation

do_plots = 1;
do_playback = 0;

%% Perform binaural rendering using the SMATBIN filters

% Recalculate NFFT_conv 
%NFFT_conv = size(drirs_simulated, 2) + size(smatbin_l, 2) - 1;
NFFT_conv =  NFFT_proc/2; 

smatbin_TF_conv_l = fft(smatbin_l, NFFT_conv, 2);
smatbin_TF_conv_l = smatbin_TF_conv_l(:, 1:round(size(smatbin_TF_conv_l, 2)/2 +1), :);

smatbin_TF_conv_r = fft(smatbin_r, NFFT_conv, 2);
smatbin_TF_conv_r = smatbin_TF_conv_r(:, 1:round(size(smatbin_TF_conv_r, 2)/2 +1), :);

DRTFs_sim_conv = fft(drirs_simulated, NFFT_conv, 2);
DRTFs_sim_conv = DRTFs_sim_conv(:, 1:round(size(DRTFs_sim_conv, 2)/2 +1));

BRTFs_smatbin_l = zeros(1, size(smatbin_TF_conv_l, 2));
BRTFs_smatbin_r = zeros(1, size(smatbin_TF_conv_r, 2));

%Convolution in frequency domain + summing
for ch_idx = 1:size(smatbin_l, 1)
    BRTFs_smatbin_l = BRTFs_smatbin_l + smatbin_TF_conv_l(ch_idx, :, head_or_idx) .* DRTFs_sim_conv(ch_idx, :);
    BRTFs_smatbin_r = BRTFs_smatbin_r + smatbin_TF_conv_r(ch_idx, :, head_or_idx) .* DRTFs_sim_conv(ch_idx, :);
end

brirs_smatbin = cat(3, sofia_tdt(BRTFs_smatbin_l), sofia_tdt(BRTFs_smatbin_r));

%% Perform binaural rendering with SOFiA (virtual loudspeaker decoding on Lebedev grid)
DRTFs_simulated = fft(drirs_simulated, NFFT_proc, 2);
DRTFs_simulated = DRTFs_simulated(:, 1:round(size(DRTFs_simulated, 2)/2)+1);

DRTFs_sma_nm = sofia_stc(N_grid, DRTFs_simulated, grid_data_sma);

rf_fullspec = [info.radial_filters, conj(info.radial_filters(:, size(info.radial_filters, 2)-1:-1:2, :))];
rf_fullspec(:, end/2+1) = real(rf_fullspec(:, size(rf_fullspec, 2)/2+1));

rf_irs = real(ifft(rf_fullspec, [], 2));
radial_filters = fft(rf_irs, NFFT_proc, 2);
radial_filters = radial_filters(:, 1:round(size(radial_filters, 2)/2)+1);

ncGap = 0;
[BRTFs_sofia_l, BRTFs_sofia_r] = sofia_binauralX(DRTFs_sma_nm, radial_filters, head_or_eval(1), 1, 1, 2048, ncGap);

brirs_sofia = cat(3, sofia_tdt(BRTFs_sofia_l), sofia_tdt(BRTFs_sofia_r));
%%
brirs_smatbin_ = brirs_smatbin;
brirs_sofia_ = brirs_sofia;
%% Plots
brirs_smatbin = brirs_smatbin_;
brirs_sofia = brirs_sofia_;
if do_plots
    smooth_plot = 0;
    NFFT_plot = 2^14;
    
    % Align BRIRs
    onset = AKonsetDetect(squeeze(brirs_sofia(1, :, 1)).', 20, -20);
    onset = ceil(mean([onset, AKonsetDetect(squeeze(brirs_sofia(1, :, 2)).', 20, -20)]))-16;
    brirs_sofia = circshift(brirs_sofia, -onset, 2);

    onset = AKonsetDetect(squeeze(brirs_smatbin(1, :, 1)).', 20, -20);
    onset = ceil(mean([onset, AKonsetDetect(squeeze(brirs_smatbin(1, :, 2)).', 20, -20)]))-16;
    brirs_smatbin = circshift(brirs_smatbin, -onset, 2);
    
    % Bring to same length
    if size(brirs_smatbin, 2) <= size(brirs_sofia, 2) 
        brirs_sofia = brirs_sofia(:, 1:size(brirs_smatbin, 2), :);
    else
        brirs_sofia = [brirs_sofia, zeros(size(brirs_sofia, 1), size(brirs_smatbin, 2)-size(brirs_sofia, 2), size(brirs_sofia, 3))];
    end
    
    % peak normalize
    brirs_sofia = brirs_sofia ./ max(max(abs(brirs_sofia)));
    brirs_smatbin = brirs_smatbin ./ max(max(abs(brirs_smatbin)));
    
    % Some plot parameters
    gap = [0.08, 0.07];
    marg_h = [0.1, 0.08]; 
    marg_w = [0.08, 0.03];

    gray = [0.7 0.7 0.7];
    blue = [12/255, 94/255, 156/255];
    red = [156/255, 27/255, 12/255];
    yellow = [249/255, 248/255, 113/255];
    cyan = [0/255, 173/255, 162/255];
    lineWidthBack = 0.5;
    lineWidthTop = 1.25;
    fontSize = 9;
    
    BRTFs_sofia = fft(brirs_sofia, NFFT_plot, 2);
    BRTFs_sofia = BRTFs_sofia(:, 1:round(size(BRTFs_sofia, 2)/2+1), :);
    
    BRTFs_smatbin = fft(brirs_smatbin, NFFT_plot, 2);
    BRTFs_smatbin = BRTFs_smatbin(:, 1:round(size(BRTFs_smatbin, 2)/2+1), :);
    
    if smooth_plot
        BRTFs_sofia_plot = AKfractOctSmooth(squeeze(BRTFs_sofia(1, :, 1)).', 'amp', fs, 6).';
        BRTFs_smatbin_plot = AKfractOctSmooth(squeeze(BRTFs_smatbin(1, :, 1)).', 'amp', fs, 6).';
    else
        BRTFs_sofia_plot = squeeze(BRTFs_sofia(1, :, 1));
        BRTFs_smatbin_plot = squeeze(BRTFs_smatbin(1, :, 1));
    end

    close all;
    figure
    subtightplot(2, 2, 1, gap, marg_h, marg_w)
        plot(abs(brirs_sofia(1, :, 1))-abs(brirs_smatbin(1, :, 1)), 'Color', gray);
        hold on;
        plot(abs(brirs_sofia(1, :, 1)), 'Color', red);
        plot(-abs(brirs_smatbin(1, :, 1)), 'Color', blue);
        title('Left')
        legend('Difference', 'SOFIA', 'SMATBIN')
        xlim([0 NFFT_pw])
        ylim([-0.5 0.5])
        xlabel('Samples')
        ylim([-1 1])
        
    subtightplot(2, 2, 2, gap, marg_h, marg_w)
        plot(abs(brirs_sofia(1, :, 2))-abs(brirs_smatbin(1, :, 2)), 'Color', gray);
        hold on;
        plot(abs(brirs_sofia(1, :, 2)), 'Color', red);
        hold on;
        plot(-abs(brirs_smatbin(1, :, 2)), 'Color', blue)
        title('Right')
        legend('Difference', 'SOFIA', 'SMATBIN')
        xlim([0 NFFT_pw])
        xlabel('Samples')
        ylim([-1 1])
    
    subtightplot(2, 2, [3, 4], gap, marg_h, marg_w)
        semilogx(linspace(5*eps, fs/2, size(BRTFs_sofia_plot, 2)), ...
                 20*log10(abs(BRTFs_sofia_plot.'))-20*log10(abs(BRTFs_smatbin_plot.')), ...
                 'LineWidth', 2, 'Color', gray)
        hold on;    
        semilogx(linspace(5*eps, fs/2, size(BRTFs_sofia_plot, 2)), ...
                 20*log10(abs(BRTFs_sofia_plot.')), 'LineWidth', 2, 'Color', red)
        
        semilogx(linspace(5*eps, fs/2, size(BRTFs_smatbin_plot, 2)), ...
                 20*log10(abs(BRTFs_smatbin_plot.')), 'LineWidth', 2, 'Color', blue, 'LineStyle', '--')
        grid on;
        xlim([20 20000])
        xticks([100 1000 10000 20000])
        xticklabels({'100', '1k', '10k', '20k'});
        legend('Difference', 'SOFIA', 'SMATBIN', 'Location', 'SouthEast')
        
        xlabel('Frequency in Hz')
        ylabel('Mmagnitude in dB')
        
    if export_plots
        % Generate path
        path = 'demo_exports';
        if ~exist(path, 'dir')
            mkdir(path);
        end
        % Export figure
        set(gcf, 'Color', 'white');
        export_fig(gcf, sprintf('%s/Simulated_pw_smatbin_%d', path, smatbin_filter_length), ...
            '-pdf','-painters','-transparent');
    end
end

%% Playback

if do_playback
    [testSignal, fs_] = audioread('MARA_LE_DRUMS.wav');
    if fs_ ~= fs
        error('Sampling rates mismatch');
    end
    bin_sound_direct = [conv(testSignal,brirs_sofia(1, :, 1)), ...
                        conv(testSignal,brirs_sofia(1, :, 2))];


    bin_sound_sma_fir = [conv(testSignal,brirs_smatbin(1, :, 1)), ...
                         conv(testSignal,brirs_smatbin(1, :, 2))];

    soundsc(bin_sound_sma_fir, fs); pause(3); clear sound;
    soundsc(bin_sound_direct, fs); pause(3); clear sound;
end
