function gen_ssr_wav(data, filename, nBits, fs)
% write a 720 channel wav file for auralization in the sound scape renderer
% in binaural synthesis mode.
%
% Parameter:
% ------------
%
% data: [M x N x Ears] M: head orientations (has to be 360)
%                      N: samples
%
% Returns:
% ------------
%
% (C) 06/2019 Tim Luebeck
% latest update: 01.02.2021
%
    if nargin < 3
        nBits = 16;
    end
    if nargin < 4
        fs = 48000;
    end
    
    ir_array = zeros( size(data, 2), 720);
    inc = 1;
                
    for node_idx = 1:360
        if ~mod(node_idx, 10)
            fprintf('|');
        end
        
        ir_array(:, inc) = data(node_idx, :, 1);
        inc = inc + 1;
        ir_array(:, inc) =  data(node_idx, :, 2);
        inc = inc + 1;
    end
    fprintf('\n\n');
    ir_array = 0.99 * ir_array ./ max(max(abs(ir_array)));
    audiowrite_multi(ir_array, fs, nBits, sprintf('%s.wav', filename));
    disp(['Sucessfully generated SSR file: ', filename])
end

