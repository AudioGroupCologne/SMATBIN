% SMATBIN: Spherical Microphone Array To Binaural 
%
% Function calc_smatbin_filter
% ----------------------------
%
% Function to calculate SMATBIN filters by successively applying a unit impulse to
% each microphone channel of the SMA, as described in [1]. In this implementation,
% the SOFiA rendering chain applying the virtual loudspeaker approach is used [2]
%
% Parameter:
% ------------
%    grid_data_sma: [Q X {Az, Col, Weights}] Q number of sampling points 
%                                            Az and Col in radiant
%                                            Weights are optional
%    N_grid:            Spatial order of the sampling scheme
%    radius:            Radius of the spherical microphone array
%    filter_length:     Final length of the SMATBIN filters
%    fs:                Sampling rate 
%    head_orientations: Head orientations for which the SMATBIN filters are
%                       calculated. [Az1, Col1; Az2, Col2; ... AzN, ColN] in radiant
%    array_config:      SOFiA-style ID specifying the array body: 
%                       0: Open sphere with pressure transducers, 
%                       2: Rigid sphere with pressure transducers {default: 2}          
%    soft_limit:        Radial filter soft-limit gain in dB {default: 20}
%    sofia_hp:          Boolean to en/disable a 30 Hz high-pass {default: 0} 
%
%
% Returns:
% -----------
%
%    smatbin_l, smatbin_r: SMATBIN filters for left and right ear
%    info: Struct with further details about the filters
%       info.radial_filters:      Radial filters applied to the SH signals
%       info.array_config:        Rigid/open sphere
%       info.radius:              Array radius
%       info.soft_limit:          Radial filter soft-limit
%       info.sofia_hp:            Boolean of 30Hz high-pass
%       info.radial_filter_delay: Latency in samples caused by linear phase radial filters
%       info.filter_proc_length:  Internal FFT size for binaural rendering 
%
%
% Dependencies:
% -------------
% SOFiA toolbox: https://github.com/AudioGroupCologne/SOFiA
%
%
% References:
% -------------
% [1] J.M. Arend, T. Lübeck, and C. Pörschmann, 
%     "Efficient binaural rendering of spherical microphone array data by linear filtering", 
%     Submitted for publication
% [2] Bernschuetz, B., Poerschmann, C., Spors, S., & Weinzierl, S. (2011),
%     "SOFiA Sound Field Analysis Toolbox. Proceedings of the International
%     Conference on Spatial Audio (ICSA)"
%
% (C) 2021 by TL & JMA, Tim Lübeck & Johannes M. Arend
%                       TH Köln - University of Applied Sciences
%                       Institute of Communications Engineering
%                       Department of Acoustics and Audio Signal Processing 
%%
function [smatbin_l, smatbin_r, info] = calc_smatbin_filter(grid_data_sma, N_grid, radius, filter_length, fs, ...
                                                                       head_orientations, array_config, soft_limit, sofia_hp)
if nargin < 6
    error('specify at least the grid data, N_grid, array radius, a desired filterlength, a sampling rate, and a set of head orientations')
end
if nargin < 7
    array_config = 2; % ac = 2 is for rigid sphere arrays with pressure transducers in SOFiA    
end
if nargin < 8
    soft_limit = 20; % Soft-limit for radial filters in dB
end
if nargin < 9
    sofia_hp = 0;
end

filter_proc_length = 8192; % Default internal filter-processing length 
if filter_length > filter_proc_length
    filter_proc_length = filter_length;
end
NFFT = filter_proc_length * 4; % Default internal FFT size
dirac_offset = 1;
nc_gap = 0;
DRTFs = zeros(size(grid_data_sma, 1), round(NFFT/2) +1);

dirac_ir.impulseResponses = zeros(1,(size(DRTFs, 2)-1)*2);
dirac_ir.impulseResponses(1+dirac_offset)=1;

% Fill SOFiA time domain struct for time-frequency Fourier transform
dirac_ir.FS = fs;
dirac_ir.averageAirTemp = 22;
dirac_ir.radius = radius;
[dirac_TF, kr] = sofia_fdt(dirac_ir, 1);

radial_filters = sofia_mf(N_grid, kr, array_config, soft_limit, 0);
[radial_filters, ~, rf_delay] = sofia_rfi(radial_filters, NFFT/4);

for ch_idx = 1:size(grid_data_sma, 1)
    fprintf('Apply dirac to microphone %d, at az: %f, col: %.2f', ch_idx, rad2deg(grid_data_sma(ch_idx,1)), rad2deg(grid_data_sma(ch_idx,2)));
    DRTFs_ch = DRTFs;
    DRTFs_ch(ch_idx, :) = dirac_TF;  
  
    % Perform binaural rendering using sofia_binauralX (virtual loudspeaker approach)
    DRTFs_nm = sofia_stc(N_grid, DRTFs_ch, grid_data_sma);
    [BRTF_l, BRTF_r] = sofia_binauralX(DRTFs_nm, radial_filters, head_orientations(:, 1), 1, sofia_hp, 2048, nc_gap);
    
    % Store filters in time domain
    smatbin_l(ch_idx, :, :) = sofia_tdt(BRTF_l).';
    smatbin_r(ch_idx, :, :) = sofia_tdt(BRTF_r).';
end

% Apply shift with respect to radial filter latency and Dirac offset
smatbin_l = circshift(smatbin_l, -(rf_delay + dirac_offset), 2);
smatbin_r = circshift(smatbin_r, -(rf_delay + dirac_offset), 2);

% Calculate window
headWinLength = 8;
tailWinLength = round(filter_length/4);
win_head = 0.5 * (1 - cos(2*pi*(0:headWinLength*2-1)'/(headWinLength*2-1)));
win_tail = 0.5 * (1 - cos(2*pi*(0:tailWinLength*2-1)'/(tailWinLength*2-1)));
winFkt = [win_head(1:headWinLength, :); ones(filter_length-headWinLength-tailWinLength, 1); win_tail(tailWinLength+1:end, :)];

% Truncate
smatbin_l = smatbin_l(:, 1:filter_length, :);
smatbin_r = smatbin_r(:, 1:filter_length, :);

% Apply window 
for k = 1:size(smatbin_l,3)
    smatbin_l(:, :, k) = smatbin_l(:, :, k) .* winFkt.';
    smatbin_r(:, :, k) = smatbin_r(:, :, k) .* winFkt.';
end 

% Save and return SMATBIN info struct
info.radial_filters = radial_filters;
info.array_config = array_config;
info.radius = radius;
info.soft_limit = soft_limit;
info.sofia_hp = sofia_hp;
info.radial_filter_delay = rf_delay; 
info.filter_proc_length = filter_proc_length;
end