function stream = simulate_sma_stream(drirs, anechoic_signal, fs, filename, multichannel, length, export, nBit)
%
% convole a sample of anechoic audio content with spherical microphone
% array impulse responeses, for offline simulation of array stream.
%
% Parameter:
% ------------
%
% data: [Q x N]  Q: number of microphones
%                N: samples
% anechoic_signal
%
% fs: sampling rate
%
% filename: filename for export
% 
% multichannel: multichannel (true), or one wave files for each SMA channel (false) 
%               {default: true}
%
% length: lenngth of audio snipped in seconds {default: 2}
% 
% export: export as wav files or just as a return value. {default: false}
%
% Returns:
% ------------
%
% (C) 06/2019 Tim Luebeck
% latest update: 01.02.2021
%
if nargin < 5
   multichannel = 0;
end
if nargin < 6
   length = 2;
end
if nargin < 7
   export = 0;
end
if nargin < 8
    nBit = 16;
end
fprintf('Simulate a spherical microphone array stream.\n')
if multichannel 
    fprintf('- as mutlichannel wav.\n'); 
end
fprintf('- length = %d seconds\n', length);
fprintf('- fs     = %d.\n', fs); 
fprintf('- nBits  = %d.\n', nBit); 
length_samples = length*fs;
anechoic_signal = anechoic_signal(:, 1:length_samples);
NFFT = size(anechoic_signal, 2) + size(drirs, 2)-1;

%% colvolve anechoic test signal with array impulse responses
DRTFs = fft(drirs, NFFT, 2);
DRTFs = DRTFs(:, 1:round(size(DRTFs, 2)/2 +1), :);

signal_TF = fft(anechoic_signal, NFFT, 2);
signal_TF = signal_TF(:, 1:round(size(signal_TF, 2)/2 +1), :);

stream_TF = DRTFs .* signal_TF;

stream_TF = [stream_TF, conj(stream_TF(:, end-1:-1:2, :))]; 
stream = ifft(stream_TF, [], 2, 'symmetric');
stream = stream ./ max(max(max(abs(stream)))) .* 0.99;

if export
   if ~multichannel
       for ch = 1:size(drirs, 1)
           audiowrite(sprintf('%s_ch%d_%dbit.wav', filename, ch, nBit), stream(ch, 1:length_samples).', fs, 'BitsPerSample', nBit);
       end
   else
        audiowrite_multi(stream(:, 1:length_samples).', fs, nBit, sprintf('%s_%dbit.wav', filename, nBit))
   end
end
end