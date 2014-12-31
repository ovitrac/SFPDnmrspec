function print_eps(resolution,filename,pathstr)
%PRINT_EPS creates an EPS file at resolution dpi  with an embedded TIFF
%   print_eps(resolution,file[,path])
%       resolution: dpi value (default=300)
%         filename: filename (full or relative)
%             path: full or relative path
%
%   See: export_fig, print2eps

% INRA\MS 2.1 - 11/06/09 - Olivier Vitrac - rev. 02/01/13

% revision
% 29/09/09 update help
% 19/12/12 change -depsc to -depsc2
% 31/12/12 use print2eps if available (its main interest is to use fix_lines)
% 02/01/13 use export_fig when available (fonts are embedded)

% arg check
resolution_default = 300;
if nargin<3, pathstr = ''; end
if ~isempty(pathstr), 
    [~,name,ext] = fileparts(filename);
else
    [pathstr,name,ext] = fileparts(filename);
end
if isempty(resolution), resolution = resolution_default; end
if isempty(pathstr), pathstr = cd; end
if ~exist(pathstr,'dir'), error('the path ''%s'' does not exist',pathstr), end
filename = [name ext];

% do print
if ispc && ~isempty(which('pdftops.exe')) && ~isempty(which('pdftops')) && ~isempty('export_fig')
    export_fig(fullfile(pathstr,filename),'-eps','-q101',sprintf('-r%d',resolution))
else
    if isempty(which('print2eps'))
        print ('-depsc2','-tiff',['-r' int2str(resolution)],fullfile(pathstr,filename))
    else
        print2eps(fullfile(pathstr,filename),gcf,'-tiff',['-r' int2str(resolution)])
    end
end