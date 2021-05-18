% SMATBIN: Spherical Microphone Array To Binaural 
%
% Script DEMO2_measured_sma_data
% ------------------------------
%
% DEMO 2: Comparison of binaural rendering of measured spherical microphone array
%         data using SMATBIN filters and the SOFiA rendering chain 
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

% Define array configuartions:
c              = 343;    
fs             = 48000;     
smatbin_filter_length = 2048;   % Length of SMATBIN filters    
soft_limit     = 20;            % Soft-limit for radial filters in dB
sofia_hp       = 1;             % Apply SOFiA's default high-pass at 30 Hz 

% Load array impulse respones and specify the array geometry
sfob = SOFAload('src/DRIR_CR1_VSA_50RS_L.sofa');
drirs = squeeze(sfob.Data.IR);
N_grid = 5;
[grid_data_sma] = get_sampling_grid('lebedev', N_grid);

radius = 0.0875;
array_body = 2;  % 2 for rigid sphear array in SOFiA      

% Head orientations 
head_orientations = get_sampling_grid('horizontal');
head_orientations = head_orientations(1:90:end, :);

%% Calculate SMATBIN filters for the SMA configuration

[smatbin_l, smatbin_r, info] = calc_smatbin_filter(grid_data_sma, N_grid, radius, ...
                                        smatbin_filter_length, fs, head_orientations, array_body, soft_limit, sofia_hp);

%% Evaluation parameters

head_or_eval = [0, pi/2]; % Head orientation for evaluation
head_or_idx = get_nn_idx(head_orientations, head_or_eval(1), head_or_eval(2)); % idx of head orientation for evaluation

do_plots = 1;
export_plots = 1;
do_playback = 0;

%% Perform binaural rendering using the SMATBIN filters

% Recalculate NFFT_conv 
NFFT_conv = size(drirs, 2) + size(smatbin_l, 2) - 1;

bin_smatbin_TF_conv_l = fft(smatbin_l, NFFT_conv, 2);
bin_smatbin_TF_conv_l = bin_smatbin_TF_conv_l(:, 1:round(size(bin_smatbin_TF_conv_l, 2)/2 +1), :);

bin_smatbin_TF_conv_r = fft(smatbin_r, NFFT_conv, 2);
bin_smatbin_TF_conv_r = bin_smatbin_TF_conv_r(:, 1:round(size(bin_smatbin_TF_conv_r, 2)/2 +1), :);

DRTFs_conv = fft(drirs, NFFT_conv, 2);
DRTFs_conv = DRTFs_conv(:, 1:round(size(DRTFs_conv, 2)/2 + 1));

BRTFs_smatbin_l = zeros(1, size(bin_smatbin_TF_conv_l, 2));
BRTFs_smatbin_r = zeros(1, size(bin_smatbin_TF_conv_l, 2));

%Convolution in frequency domain + summing
for ch_idx = 1:size(smatbin_l, 1)
    BRTFs_smatbin_l = BRTFs_smatbin_l + bin_smatbin_TF_conv_l(ch_idx, :, head_or_idx) .* DRTFs_conv(ch_idx, :);
    BRTFs_smatbin_r = BRTFs_smatbin_r + bin_smatbin_TF_conv_r(ch_idx, :, head_or_idx) .* DRTFs_conv(ch_idx, :);
end
brirs_smatbin_l = sofia_tdt(BRTFs_smatbin_l);  
brirs_smatbin_r = sofia_tdt(BRTFs_smatbin_r);

brirs_smatbin = cat(3, sofia_tdt(BRTFs_smatbin_l), sofia_tdt(BRTFs_smatbin_r));

%% Perform binaural rendering with SOFiA (virtual loudspeaker decoding on Lebedev grid)

% Recalculate NFFT_proc
NFFT_proc = 2*size(drirs, 2) + 2*smatbin_filter_length; 
DRTFs = fft(drirs, NFFT_proc, 2);
DRTFs = DRTFs(:, 1:round(end/2)+1);

DRTFs_sma_nm = sofia_stc(N_grid, DRTFs, grid_data_sma);

rf_fullspec = [info.radial_filters, conj(info.radial_filters(:, size(info.radial_filters, 2)-1:-1:2, :))];
rf_fullspec(:, end/2+1) = real(rf_fullspec(:, size(rf_fullspec, 2)/2+1));

rf_irs = real(ifft(rf_fullspec, [], 2));
radial_filters = fft(rf_irs, NFFT_proc, 2);
radial_filters = radial_filters(:, 1:round(size(radial_filters, 2)/2)+1);

[BRTFs_sofia_l, BRTFs_sofia_r] = sofia_binauralX(DRTFs_sma_nm, radial_filters, head_or_eval(1), 1, sofia_hp); %Lebedev composite grid (virtual loudspeakers)

%% Plots

% Bring to same length
brirs_sofia = cat(3, sofia_tdt(BRTFs_sofia_l), sofia_tdt(BRTFs_sofia_r));
if size(brirs_smatbin, 2) <= size(brirs_sofia, 2) 
    brirs_sofia = brirs_sofia(:, 1:size(brirs_smatbin, 2), :);
else
    brirs_sofia = [brirs_sofia, zeros(size(brirs_sofia, 1), size(brirs_smatbin, 2)-size(brirs_sofia, 2), size(brirs_sofia, 3))];
end

% Align BRIRs
onset = AKonsetDetect(squeeze(brirs_sofia(1, :, 1)).', 20, -20);
onset = ceil(mean([onset, AKonsetDetect(squeeze(brirs_sofia(1, :, 2)).', 20, -20)]))-16;
brirs_sofia = circshift(brirs_sofia, -onset, 2);

onset = AKonsetDetect(squeeze(brirs_smatbin(1, :, 1)).', 20, -20);
onset = ceil(mean([onset, AKonsetDetect(squeeze(brirs_smatbin(1, :, 2)).', 20, -20)]))-16;
brirs_smatbin = circshift(brirs_smatbin, -onset, 2);

smooth_plot = 1;
if do_plots
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
    
    BRTFs_sofia = fft(brirs_sofia, [], 2);
    BRTFs_sofia = BRTFs_sofia(:, 1:round(size(BRTFs_sofia, 2)/2+1), :);
    BRTFs_smatbin = fft(brirs_smatbin, [], 2);
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
        xlim([0 3000])
        ylim([-0.5 0.5])
        xlabel('Samples')
        
    subtightplot(2, 2, 2, gap, marg_h, marg_w)
        plot(abs(brirs_sofia(1, :, 2))-abs(brirs_smatbin(1, :, 2)), 'Color', gray);
        hold on;
        plot(abs(brirs_sofia(1, :, 2)), 'Color', red);
        hold on;
        plot(-abs(brirs_smatbin(1, :, 2)), 'Color', blue)
        title('Right')
        legend('Difference', 'SOFIA', 'SMATBIN')
        xlim([0 3000])
        xlabel('Samples')
    
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
        export_fig(gcf, sprintf('%s/Measured_sma_smatbin_%d', path, smatbin_filter_length), ...
            '-pdf','-painters','-transparent');
    end
end

%% Playback

if do_playback
    [testSignal, fs_] = audioread('MARA_LE_DRUMS.wav');
    if fs_ ~= fs
        error('Sampling rates mismatch');
    end
    bin_sound_direct = [conv(testSignal, brirs_sofia(1, :, 1)), ...
                        conv(testSignal, brirs_sofia(1, :, 2))];

       
    bin_sound_smatbin = [conv(testSignal, brirs_smatbin(1, :, 1)), ...
                         conv(testSignal, brirs_smatbin(1, :, 2))];
       
    soundsc(bin_sound_smatbin, fs); pause(3); clear sound;
    soundsc(bin_sound_direct, fs); pause(3); clear sound;
end
