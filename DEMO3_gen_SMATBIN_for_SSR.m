% SMATBIN: Spherical Microphone Array To Binaural 
%
% Script DEMO3_gen_SMATBIN_for_SSR
% --------------------------------
%
% DEMO 3: SMATBIN decoding using the SoundScape Renderer: Calculate SMATBIN 
% filters for a 19 channel Zylia ZM-1 SMA and export for dynamic binaural 
% synthesis using the SoundScape Renderer. Optionally, generate a
% simulated spherical microphone array stream for offline testing
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
%
%%
clear all; close all; clc

gen_fake_stream = 1; % If true, a simulated SMA stream is generated for offline testing 
                     
%% Define array configuartions

c              = 343;    
fs             = 48000;     
smatbin_filter_length = 2048;   % Length of SMATBIN filters    
soft_limit     = 20;            % Soft-limit for radial filters in dB
sofia_hp       = 1;             % Apply SOFiA's default high-pass at 30 Hz 

% Specify array geometry
array_body = 2;  % 2 for rigid sphear array in SOFiA  

% Get Zylia ZM-1 sampling grid, N, and radius
[grid_data_sma, ~, N_grid, radius] = get_sampling_grid('zylia');
%[grid_data_sma, ~, N_grid, radius] = get_sampling_grid('extremal', 12);
%radius = 0.0875;

% Define head orientations to be rendered
head_orientations = get_sampling_grid('horizontal'); % Render 360 horizontal head orientations

%% Calculate SMATBIN filters for the SMA configuration

[smatbin_l, smatbin_r, info] = calc_smatbin_filter(grid_data_sma, N_grid, radius, ...
                                        smatbin_filter_length, fs, head_orientations, array_body, soft_limit, sofia_hp);

%% SSR export 

% Generate path
path = 'demo_exports';
if ~exist(path, 'dir')
    mkdir(path);
end

% Export as SSR filters:
smatbin_ssr = cat(4, smatbin_l, smatbin_r); % Bring to shape [Mic x samples x head orientation x ears]
export_sma_firs_SSR(smatbin_ssr, fs, sprintf('%s/demo_N%d', path, N_grid));

%% Simulate a sound field with plane waves and convolve with test signal
if gen_fake_stream
    NFFTpw = 1024;
    az   = 0*(pi/180);    
    col  = 90*(pi/180);   
    ds = NFFTpw/2/fs;          

    [DRTFs_simulated, ~] = sofia_swg(radius, grid_data_sma, array_body, fs, ...
                             NFFTpw, az, col, 100, ds);    

    DRTFs_simulated = [DRTFs_simulated, conj(DRTFs_simulated(:, end-1:-1:2, :))];
    DRTFs_simulated(:, NFFTpw/2+1) = real(DRTFs_simulated(:, NFFTpw/2+1));
    drirs_simulated = real(ifft(DRTFs_simulated, [], 2));

    % Generate an array stream for SSR evaluation
    [anechoic_signal, ~] = audioread('MARA_LE_DRUMS.wav');
    simulate_sma_stream(drirs_simulated, anechoic_signal.', fs, ...
        sprintf('%s/sma_stream_drums_N%d', path, N_grid), 1, 3, 1, 16);
end






