function lags = nmrfindlags(y1,y2,varargin)
%NMRFINDLAGS find lags corresponding to maxima of correlation when 2 signals are correlated
% synatx lags = nmrfindlags(y1,y2,'maxlagsx',[value]...)
%
% INPUTS
%            y1: m1x1 vector of y valuee (y1 is movable relative to y2)
%            y2: m2x1 vector of y values (y2 is fixed)
%            dx: step in x (default = 1)
%          nmax: maximum number of lags found (default = inf)
%      maxlagsx: acceptable maxlags in x unit 
%       maxlags: acceptable maxlags in index (default = m = max(m1,m2))
%          zero: zero value for monotonepeak (default = sqrt(eps)/1e3)
%         mfilt: filtering width for filt data for monotonepeak (default m/20)
%
% OUTPUTS
%         lags : array,(nmax x 1 as size)
%
% RMNSPEC v 0.1 - 06/05/2013 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 23/5/2013
% history
% 23/5/13 transpose xlags
%         add zero value and mfilt for monotonepeak
%         fix empty p when no peaks is found in xc -> lags corresponds to max(xc)

% default
default = struct('maxlagsx',[],'maxlags',[],'dx',[],'nmax',[],'zero',[],'mfilt',[]);

% argcheck
o = argcheck(varargin,default);
if nargin<2, error('2 arguments are required'), end
y1 = y1(:); y2 = y2(:); 
m1 = length(y1); m2 = length(y2); m = max(m1,m2);
if isempty(o.dx), o.dx = 1; end
if ~isempty(o.maxlagsx), o.maxlags = ceil(o.maxlagsx/o.dx); end
if isempty(o.maxlags), o.maxlags = m; end
if isempty(o.nmax), o.nmax = inf; end
if isempty(o.zero), o.zero = sqrt(eps)/1e3; end
if isempty(o.mfilt), o.mfilt = m/20; end

% zeropadding on the right hand side of initial signal       
y1 = [y1;zeros(m-m1,1)]; % 
y2 = [y2;zeros(m-m2,1)];      

% main
[xc,xlags] = xcorr(y1/max(y1),y2/max(y2),o.maxlags,'unbiased'); 

% search maxima of cglobal
p = monotonepeak('x',xlags','y',xc,'zero',o.zero,'mfilt',o.mfilt,'sort','descend','sortby','absheight','array'); % 'mfilt',1/50*length(xlags)
if ~isempty(p)
    nlag = min([o.nmax size(p,1)]); 
    lags = cat(1,p(1:nlag).center);
else [~,imax] = max(xc);
    lags = xlags(imax);
end
