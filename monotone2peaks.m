function peaks=monotone2peaks(I,varargin)
%MONOTONE2PEAKS find peaks using monotone with '+' and '-' arguments
% SYNTAX
%   peaks = monotone2peaks(I [,keyword1,value1,keyword2,value2,...])
%   peaks = monotone2peaks(...)
% INPUTS
%     I: mx1 array (uniform sampling)
%     pair keyword/value (default value)
%       'mfilt': 1xn filtering values (to be used with filtzero) (default = 1), 1 means no filtering
%         'ref': index of reference filetring width (default=1) (see ldist)
%           'm': resolution for base commun (X, Y) (default=1e4) (see ldist)
%      'method': method for interpolation 'interp1)(default='cubic') (see ldist)
% OUTPUTS (n=1, i.e. mfilt is a scalar)
%     peaks.If: npx1 array of filtered data according to mfilt
%     peaks.p: npx1 array of positions/indices of peaks ranked from the tallest to the lowest, where np is the number of peaks
%     peaks.a: npx1 array of corresponding heights (max of values returned by mononotone with '+' and '-' arguments
%     peaks.rr: relative rank of the peak (i/np for the peak i)
%     peaks.pstart: position of the beginning of peaks 
%     peaks.pstop: position of the end of peaks 
%     peaks.preject: percentage of excluded peaks        
%     peaks.areject: amplitude at rejection percentage 

% OUTPUTS (n>1, i.e. mfilt is an array)
%     peaks(i).If: mxn array of filtered data of I according to mfilt
%     peaks(i).p: nx1 cell array where p{i} is a npix1 array of positions/indices
%     peaks(i).a: nx1 cell array where a{i} is a npix1 array of heights
%     peaks(i).rr: nx1 cell array where rr{i} is a npix1 array of relative ranks
%     peaks(i).pstart: position of the beginning of peaks 
%     peaks(i).pstop: position of the end of peaks 
%     peaks(i).preject: nx1 percentage of excluded peaks        
%     peaks(i).areject: nx1 amplitude at rejection percentage
    
%
% See also: filtzero, monotone, ldist
%
% Example
%   if isempty(find_path_toolbox('rmnspec')), error('install first the toolbox rmnspec'), end
%   [dbpur,dbxpur] = nmrloadascii('path', fullfile(find_path_toolbox('rmnspec'),'data_pur'));
%   mfilt = [1 10 20 30 50];
%   peaks=monotone2peaks(dbxpur.I(:,1),'mfilt',mfilt);
%   figure
%   hp = plotpub({peaks.rr},{peaks.a},'marker','none','linewidth',1,'linestyle','-','color',cbrewer('div','RdYlGn',length(peaks)));
%   legendpub(hp,arrayfun(@(m) sprintf('m_{filt}=%d',m),mfilt,'UniformOutput',false),1);
%   xlabel('relative rank/probability'), ylabel('peak height')
%   figure
%   i=4
%   plot(dbxpur.ppm,peaks(i).If)
%   arrayfun(@(p,rr) text(dbxpur.ppm(p),peaks(i).If(p),sprintf('%0.2g',rr),'fontsize',8,'rotation',45),peaks(i).p(1:100),peaks(i).rr(1:100))
%   
% Example (debug 05/10/2012) good S/N ratio
%   if isempty(find_path_toolbox('rmnspec')), error('install first the toolbox rmnspec'), end
%   [dbpur,dbxpur] = nmrloadascii('path', fullfile(find_path_toolbox('rmnspec'),'data_pur'));
%   mfilt = 1:3:60;
%   peaks=monotone2peaks(dbpur.Stearicacid.I,'mfilt',mfilt);
%   figure
%   hp = plotpub({peaks.rr},{peaks.a},'marker','none','linewidth',1,'linestyle','-','color',cbrewer('div','RdYlGn',length(peaks)));
%   legendpub(hp,arrayfun(@(m) sprintf('m_{filt}=%d',m),mfilt,'UniformOutput',false),1);
%   xlabel('relative rank/probability'), ylabel('peak height')
%
%Example: good signal to noise ratio
%   if isempty(find_path_toolbox('rmnspec')), error('install first the toolbox rmnspec'), end
%   [dbpur,dbxpur] = nmrloadascii('path', fullfile(find_path_toolbox('rmnspec'),'data_pur'));
%   mfilt = 1:3:60;
%   peaks=monotone2peaks(dbpur.Stearicacid400.I,'mfilt',mfilt);
%   figure
%   hp = plotpub({peaks.rr},{peaks.a},'marker','none','linewidth',1,'linestyle','-','color',cbrewer('div','RdYlGn',length(peaks)));
%   legendpub(hp,arrayfun(@(m) sprintf('m_{filt}=%d',m),mfilt,'UniformOutput',false),1);
%   xlabel('relative rank/probability'), ylabel('peak height')
 % 
% Example: bad signal to noise ratio
%   if isempty(find_path_toolbox('rmnspec')), error('install first the toolbox rmnspec'), end
%   [dbpur,dbxpur] = nmrloadascii('path', fullfile(find_path_toolbox('rmnspec'),'data_pur'));
%   mfilt = 1:3:60;
%   peaks=monotone2peaks(dbpur.Stearicacid100.I,'mfilt',mfilt);
%   figure
%   hp = plotpub({peaks.rr},{peaks.a},'marker','none','linewidth',1,'linestyle','-','color',cbrewer('div','RdYlGn',length(peaks)));
%   legendpub(hp,arrayfun(@(m) sprintf('m_{filt}=%d',m),mfilt,'UniformOutput',false),1);
%   xlabel('relative rank/probability'), ylabel('peak height')
%
% Example: bad signal to noise ratio
%   if isempty(find_path_toolbox('rmnspec')), error('install first the toolbox rmnspec'), end
%   [dbpur,dbxpur] = nmrloadascii('path', fullfile(find_path_toolbox('rmnspec'),'data_pur'));
%   mfilt = 1:3:60;
%   peaks=monotone2peaks(dbpur.MBOCA1250.I,'mfilt',mfilt);
%   figure
%   hp = plotpub({peaks.rr},{peaks.a},'marker','none','linewidth',1,'linestyle','-','color',cbrewer('div','RdYlGn',length(peaks)));
%   legendpub(hp,arrayfun(@(m) sprintf('m_{filt}=%d',m),mfilt,'UniformOutput',false),1);
%   xlabel('relative rank/probability'), ylabel('peak height')
%
% RMNSPEC v 0.1 - 26/09/12 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 06/10/12
% 
% History
% 26/09/2012 help fixes
% 28/09/2012 output as a structure
% 29/09/2012 calculate preject, areject
% 01/10/2012 fix help
% 05/10/2012 add example to debug baseline problem
% 06/10/2012 add examples (different cases of signal to noise ratio)
% 16/10/2012 fix n=1, add field n
%
% default
default = struct('mfilt',[1 2 10],'ref',1,'m',1e4,'method','cubic');

% argcheck
if nargin<1, error('one argument is required'); end
if size(I,2)>1, error('I must be a column vector'), end
o = argcheck(varargin,default);
n = numel(o.mfilt);

% main
peaks = repmat(struct('mfilt',[],'n',[],'p',[],'a',[],'rr',[],'If',[],'pstart',[],'pstop',[],'preject',[],'areject',[]),n,1);
for i=1:n
    peaks(i).mfilt = o.mfilt(i);
    % filtering
    if o.mfilt(i)==1, peaks(i).If = I;
    else peaks(i).If = filtzero(I,o.mfilt(i));
    end
    % search monotone sections
    [pp,lp,ap] = monotone(peaks(i).If,'+');
    [pm,lm,am] = monotone(peaks(i).If,'-');
    pp = pp+lp-1;
    % populate peaks (from left or from right)
    p0 = intersect(pp,pm);  % position of peaks (union to be avoided, to prevent NaN at start and stop)
    n0 = length(p0);    % number of peaks
    a0 = NaN(n0,2);     % amplitude of peaks
    [start,stop] = deal(NaN(n0,1));
    [~,ip,i0] = intersect(pp,p0);
    a0(i0,1) = ap(ip);
    start(i0) = pp(ip)-lp(ip)+1;
    [~,im,i0] = intersect(pm,p0);
    a0(i0,2) = -am(im);
    stop(i0) = pm(im)+lm(im)-1;
    % sort peaks by size
    [a0,ia] = sort(max(a0,[],2),'descend');
    % store results
    peaks(i).p = p0(ia);
    peaks(i).a = a0;
    peaks(i).rr = 100*(1:n0)'/n0;
    peaks(i).pstart = start(ia);
    peaks(i).pstop = stop(ia);
    peaks(i).n = length(ia);
end

% add ldist data
if n>1
    [ptmp,atmp] = ldist({peaks.rr},{peaks.a},o);
    for i=1:n
        peaks(i).preject = ptmp(i);
        peaks(i).areject = atmp(i);
    end
end
