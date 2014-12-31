function print_png(resolution,filename,filepath,options,flipon,margin,negativeon,horizontalon)
%PRINT_PNG  print active window as a PNG image
%   print_png(resolution,filename[,filepath,options])
%   print_png(resolution,filename,options)
%   print_png(resolution,filename,options [,flipon,margin,negativeon])
%       options must be followed by '-' (ex. '-opengl')
%       imfile,flipon,margin,negativeon are parameters of pngtruncateim

% INRA\MS 2.0 - 28/01/01 - Olivier Vitrac - 10/05/14

% revision
% 25/07/07 add options and fix filepath option
% 13/08/07 fix path ambiguity when both filename and filepath contain a path
% 09/09/12 add pngtruncateim
% 10/05/14 force flipon when paperorientation is set to landscape

% Default
flipon_default = false;
margin_default = 0;
negativeon_default = false;
horizontalon_default = false;

% arg check
if nargin<2, error('Two arguments are required'); end
if nargin<3, filepath = ''; end
if nargin<4, options = ''; end
cropon = (nargin>=5);
if nargin<5, flipon = []; end
if nargin<6, margin = []; end
if nargin<7, negativeon= []; end
if nargin<8, horizontalon = []; end
if isempty(flipon), flipon = flipon_default; end
if ~flipon && strcmpi(get(gcf,'PaperOrientation'),'landscape'), flipon = true; end
if isempty(margin), margin = margin_default; end
if isempty(negativeon), negativeon = negativeon_default; end
if isempty(horizontalon), horizontalon = horizontalon_default; end

% do the print
if any(filepath), 
    if filepath(1)=='-', options = filepath; filepath = ''; end
    [~,name,ext] = fileparts(filename);
else
    [filepath,name,ext] = fileparts(filename);
end
if isempty(ext), ext = '.png'; end
if isempty(filepath), filepath = cd; end
if ~exist(filepath,'dir'), error('the path ''%s'' does not exist',filepath), end
filename = [name ext];
print (['-r' int2str(resolution)],'-dpng',options, fullfile(filepath,filename))
if cropon
    pngtruncateim(fullfile(filepath,filename),flipon,margin,negativeon,horizontalon) % varargin{:})
else
    dispf('to truncate the generated PNG image, use the following code:\npngtruncateim(''%s'',0,0,0)',fullfile(filepath,filename))
end