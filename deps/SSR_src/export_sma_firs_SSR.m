function export_sma_firs_SSR(IRs, fs, name)
% export sma brirs for convolution in sound scape renderer
%
% Parameter:
% ------------
%    IRs:     filters to export for ssr in shape [Mic x samples x head orientation x ears]
%    fs:      sampling rate
%    name:    filename
%
% Returns:
% -----------
%
% Dependencies:
% -------------
% gen_ssr_asd, gen_ssr_wav
%
% 26.02.2021 Tim Luebeck

    filterset_filename = sprintf('%s_SSR_filters', name);
    if ~exist(filterset_filename)
        mkdir(filterset_filename);
    end
    fir_names = [];
    for mic = 1:size(IRs, 1)
        fir_name = sprintf('%s/smatbin_ch%d', filterset_filename, mic);
        filters = permute(squeeze(IRs(mic, :, :, :)), [2, 1, 3]);
        gen_ssr_wav(filters, fir_name, 24, fs);
        fir_names{mic} = sprintf('smatbin_ch%d.wav', mic);
    end
    gen_ssr_asd(sprintf('%s/smatbin_scene', filterset_filename), fir_names);
end


