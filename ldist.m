function [x,aout] = ldist(X,Y,varargin)
%LDIST calculates the x-position of a L corner using the cross algorithm
% SYNTAX
%   x=ldist(X,Y,[,'property',value,...)
% Property/value
%       'ref': index of reference filetring width (default=1)
%         'm': resolution for commun base  (X, Y) (default=1e4)
%    'method': method for interpolation (interp1)(default='cubic')
% 
% OUTPUT
%           x: percentage of excluded peaks
%           a: amplitude at rejection percentage
%
% Example
%   if isempty(find_path_toolbox('rmnspec')), error('install first the toolbox rmnspec'), end
%   [dbpur,dbxpur] = nmrloadascii('path', fullfile(find_path_toolbox('rmnspec'),'data_pur'),'ppmstandards',[-0.05 0.15; 1.50 1.65; 7.2 7.32;],'ppmmin',-0.5,'ppmmax',12);
%   peaks=monotone2peaks(dbxpur.I(:,1),'mfilt',1:30);
%   [preject,areject] = ldist({peaks.rr},{peaks.a}); figure, plot(preject), figure, plot(areject)

% RMNSPEC v 0.1 - 27/09/12 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 01/10/12
% 
% History
% 28/09/12 RC
% 29/09/12 fix nargin, add output a
% 01/10/12 fix help
% default
default = struct('ref',1,'m',1e4,'method','cubic');

% arg check
if nargin<2, error('2 arguments are required'), end
o = argcheck(varargin,default);
if ~iscell(X)
    if size(X,2)==1, X={X}; else X = num2cell(X,2); end
end
if ~iscell(Y)
    if size(Y,2)==1, Y={Y}; else X = num2cell(Y,2); end
end
n = length(Y);
if n<=1, dispf('WARNING: the cross algorithm requires at least two L curves'), x=NaN; return, end
if length(X)~=length(Y), X = X(ones(1,n)); end

% search bounds
xmax = cellfun(@max,X);
ymax = cellfun(@max,Y);

% main
[x,a] = deal(NaN(n,1));
for i=setdiff(1:n,o.ref)
    
    % check size compatibility
    if ~matcmp(size(X{i}),size(Y{i})), error('the sizes of X{%d} and Y{%d} do not match',i,i), end
    
    % common scale
    xbase = linspace(0,min(xmax([o.ref i])),o.m)';
    ybase = linspace(0,min(ymax([o.ref i])),o.m)';
    
    % normalized vertical distance
    Ydist = abs(interp1(X{i},Y{i},xbase,o.method)-interp1(X{o.ref},Y{o.ref},xbase,o.method))/mean(ymax([o.ref i]));
    
    % horizontal distance
    xleft = interpleft(Y{i},X{i},ybase);
    xleftref = interpleft(Y{o.ref},X{o.ref},ybase);
    [xcenter,ixcenter] = sort(mean([xleft xleftref],2),'ascend');
    Xdist = interpleft(xcenter,abs(xleft(ixcenter)-xleftref(ixcenter))/mean(xmax([o.ref i])),xbase);
    
    % search to equate vertical and horizontal distances
    d = @(x) interp1(xbase,Ydist,x,o.method)-interp1(xbase,Xdist,x,o.method); % criterion
    msearchmax = round(o.m/2); % search in the left region only
    dsearch = d(xbase(1:msearchmax));
    i0 = find((dsearch(1:end-1).*dsearch(2:end))<0,1,'last'); % last index with change of sign
    x(i) = fsolve(d,xbase(i0),optimset('Tolfun',1e-9,'TolX',1e-9,'Algorithm','trust-region-reflective','MaxIter',1e3,'Display','off')); % refinement
    a(i) = interp1(X{i},Y{i},x(i),o.method);
end

% output
if nargout>1, aout = a; end