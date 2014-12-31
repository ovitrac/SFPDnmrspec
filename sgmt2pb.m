function w = sgmt2pb(x,sgmt,varargin)
%SGMT2PB converts segments [x1min x1max;x2min x2max;x3min x3max;...] into weights using pb
% SYNTAX
%   w = sgmt2pb(x,sgmt [,keyword1,value1,keyword2,value2,...])
% INPUTS
%     x: mx1 array (scale) 
%  sgmt: nx2 array of segments coded as [x1min x1max;x2min x2max;x3min x3max;...]
%     pair keyword/value (default value)
%       'buffer' buffer width (default = 1)
%       'mode' buffer type (default = 'p')
% OUTPUT
%     w: mx1 array of weights (with value between 0 and 1)
%
% See also: pb
%
% Example
%   if isempty(find_path_toolbox('rmnspec')), error('install first the toolbox rmnspec'), end
%   [dbpur,dbxpur] = nmrloadascii('path', fullfile(find_path_toolbox('rmnspec'),'data_pur'),'ppmstandards',[-0.05 0.15; 1.50 1.65; 7.2 7.32;],'ppmmin',-0.5,'ppmmax',12);
%   sgmt = [-0.05 0.15; 1.50 1.65; 7.2 7.32];
%   w = sgmt2pb(dbxpur.ppm,sgmt,'buffer',0.05);
%   figure, plot(dbxpur.ppm,w)
%   figure, plot(dbxpur.ppm,dbxpur.I(:,1).*w)

% RMNSPEC v 0.1 - 26/09/12 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev.
% 
% History

% default
default = struct('buffer',1,'mode','p');

% argcheck
if nargin<2, error('2 arguments are required'), end
if size(x,2)>1, error('x must be a column vector'), end
o = argcheck(varargin,default);
if size(sgmt,2)~=2, error('sgmt must be a nx2 array'), end
m = size(x,1);
n = size(sgmt,1);

% main
p = sgmt'; p = sort(p(:),'ascend');
w = 1 - pb(x,p,o.buffer/4,o.mode);