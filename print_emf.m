function print_emf(resolution,fichier,chemin,options)
%PRINT_JPG  print active window as a EMF image
%   print_emf(resolution,fichier[,chemin,options])
%   print_emf(resolution,fichier,options)
%       options must be followed with '-' (ex. '-opengl')

% MS 2.1 - 28/01/01 - Olivier Vitrac - 11/04/12

% revision
% 11/04/12 release candidate derived from pint_png

if nargin<3, chemin = ''; end
if nargin<4, options = ''; end
if any(chemin), 
    if chemin(1)=='-', options = chemin; chemin = ''; end
    [~,name,ext] = fileparts(fichier);
else
    [chemin,name,ext] = fileparts(fichier);
end
if isempty(chemin), chemin = cd; end
if isempty(ext), ext = '.emf'; end
if ~exist(chemin,'dir'), error('the path ''%s'' does not exist',chemin), end
fichier = [name ext];
print (['-r' int2str(resolution)],'-dmeta','-painters',options, fullfile(chemin,fichier))