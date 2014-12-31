function out = pngtruncateim(imfile,flipon,margin,negativeon,horizontalon)
%PNGTRUNCATEIM  crop and rotate PNG already saved images
%   Syntax: pngtruncateim(imfile [,flipon,margin,negative])
%       default values: flipon = false
%                       margin = 100
%                       negativeon = false

% Migration 2.0 - 24/05/11 - INRA\Olivier Vitrac - rev. 29/03/12

% Revision history
% 27/05/11 guess extension
% 07/11/11 add negative
% 29/01/12 accept image data
% 29/03/12 fix negativeon
% 02/10/12 add horizontalon

% default values
flipon_default = false;
margin_default = 100;
negativeon_default = false;
horizontalon_default = false;

% arg check
if nargin<2, flipon = []; end
if nargin<3, margin = []; end
if nargin<4, negativeon= []; end
if nargin<5, horizontalon = []; end
if isempty(flipon), flipon = flipon_default; end
if isempty(margin), margin = margin_default; end
if isempty(negativeon), negativeon = negativeon_default; end
if isempty(horizontalon), horizontalon = horizontalon_default; end

% main section
if ischar(imfile)
    if ~exist(imfile,'file'), error('the file ''%s'' does not exist',imfile), end
    im = imread(imfile);
    nofile = false;
elseif isnumeric(imfile)
    im = imfile;
    nofile = true;
else
    error('the argument must be a valid image filename or image data')
end

siz = size(im);
if negativeon, im = 255-im; end
imb = min(im,[],3);
lim = zeros(2,2);
dimlist = 1:2;
for dim = dimlist
    lim(dim,1) = find(min(imb,[],dimlist(mod(dim,2)+1))<255,1,'first');
    lim(dim,2) = find(min(imb,[],dimlist(mod(dim,2)+1))<255,1,'last');
end
lim(:,1) = max(lim(:,1)-margin,1);
lim(:,2) = min(lim(:,2)+margin,siz(1:2)');
im = im(lim(1,1):lim(1,2),lim(2,1):lim(2,2),:);
if flipon, im = flipdim(permute(im,[2 1 3]),2); end
if negativeon, im = 255-im; end
if horizontalon && size(im,1)>size(im,2), im = flipdim(permute(im,[2 1 3]),2); end

% save final result
if nofile
    out = im;
else
    ext = uncell(regexp(imfile,'\.([^\.]+)$','tokens'));
    if isempty(ext)
        imwrite(im,imfile,'png');
    else
        imwrite(im,imfile,ext{1});
    end
end