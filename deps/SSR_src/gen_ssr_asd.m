function gen_ssr_asd(filename, BRIR_names, scene_title, hornhack, az_range, el_range)
% write a asd file for sound scape renderer
%
% Parameter:
% ------------
%
% filname: asd-filename
% BRIR_name: relative path (relative to asd filename) and file name of the BRIR wavs
%
% Returns:
% ------------
%
% (C) 06/2019 Tim Luebeck
% latest update: 01.02.2021
%
if nargin < 3
    scene_title = [];
end
if nargin < 4
    hornhack = 0;
end
if nargin < 5 && hornhack == 1
    error('Pass a azimuth range when using sofa files!')
end
if nargin < 6 && hornhack == 1
    error('Pass a elevation range when using sofa files!')
end
file = fopen(sprintf('%s.asd', filename), 'w');
fprintf(file, '<?xml version="1.0" encoding="utf-8"?>\n');
fprintf(file, '<asdf>\n');
fprintf(file, '  <header>\n');
fprintf(file, '     <name>%s</name>\n', scene_title);
fprintf(file, '     <description>\n');
fprintf(file, '     </description>\n');
fprintf(file, '  </header>\n');
fprintf(file, '  <scene_setup>\n');
fprintf(file, '    <volume>-6</volume>\n');

% calc gui positions
% SSR gui: first dimension positive to the right (Matlab convention: -y)
%          second dimension positive to the top (Matlab convention: +x)
angle_per_src = 2*pi/ length(BRIR_names);
r_gui = 3;

cur_angle = 0;
[x, y, ~] = sph2cart(cur_angle, 0, r_gui);
gui_pos = [y, -x];
for src_idx = 1:length(BRIR_names)
    fprintf(file, '    <source name="%d. %s" properties_file="%s" volume="0" mute="true">\n', ...
                    src_idx, string(BRIR_names{src_idx}), string(BRIR_names{src_idx}));
    if hornhack 
        fprintf(file, '        <properties \n', gui_pos(1), gui_pos(2));
        fprintf(file, '          brir_file="%s"\n', string(BRIR_names{src_idx}));
        fprintf(file, '          azimuth_range="%d %d"\n', az_range(1), az_range(2));
        fprintf(file, '          elevation_range="%d %d"/>\n', el_range(1), el_range(2));
    end
    
    fprintf(file, '        <position x="%f" y="%f" fixed="true"/>\n', gui_pos(1), gui_pos(2));
    fprintf(file, '    </source>\n');
    
    % update gui position
    cur_angle = cur_angle + angle_per_src;
    [x, y, ~] = sph2cart(cur_angle, 0, r_gui);
    gui_pos = [y, -x];
end
fprintf(file, '  </scene_setup>\n');
fprintf(file, '</asdf>\n');
fclose(file);
end

