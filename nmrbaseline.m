function [isbase,Ib] = nmrbaseline(I,varargin)
%NMRBASELINE fits baseline of NMR spectra
%   Syntax : [isbase,Ib] = nmrbaseline(I[,'property',value,...)
%          I : nx1 vector 
%   Property/value
%          'mfilt' : (default = 10)
%          'prejection' : percentage of excluded peaks (default: 5=5%)
%
% Output :
%          isbase : nx1 logical vector (true if I is associated to baseline)
%          Ib : nx1 vector as I with baseline removed
%
% See also: nrlmoadascii, nmrcorr, nmrsubdb
%
% Example
%   if isempty(find_path_toolbox('rmnspec')), error('install first the toolbox rmnspec'), end
%   [dbpur,dbxpur] = nmrloadascii('path', fullfile(find_path_toolbox('rmnspec'),'data_pur'),'ppmstandards',[-0.05 0.15; 1.50 1.65; 7.2 7.32;],'ppmmin',-0.5,'ppmmax',12);
%   peaks=monotone2peaks(dbxpur.I(:,1),'mfilt',1:2);
%   [isbase,Ib] = nmrbaseline(dbxpur.I(:,1),'mfilt',20,'prejection',peaks(2).preject);
%   figure,plot(dbxpur.ppm,dbxpur.I(:,1),'b-'), hold on,  plotequi(dbxpur.ppm(~isbase),dbxpur.I(~isbase,1),'r-')
%
% Example (debug 05/10/2012) good S/N ratio
%   if isempty(find_path_toolbox('rmnspec')), error('install first the toolbox rmnspec'), end
%   [dbpur,dbxpur] = nmrloadascii('path', fullfile(find_path_toolbox('rmnspec'),'data_pur'),'noprefetch');
%   peaks=monotone2peaks(dbpur.Stearicacid.I,'mfilt',[1 10]);
%   [isbase,Ib] = nmrbaseline(dbpur.Stearicacid.I,'mfilt',10,'prejection',3);
%   figure,plot(dbpur.Stearicacid.ppm,dbpur.Stearicacid.I,'b-'), hold on,  plotequi(dbpur.Stearicacid.ppm(~isbase),dbpur.Stearicacid.I(~isbase,1),'r-'), plot(dbpur.Stearicacid.ppm,peaks(2).If,'g')
%
% Example: good signal to noise ratio
%   if isempty(find_path_toolbox('rmnspec')), error('install first the toolbox rmnspec'), end
%   [dbpur,dbxpur] = nmrloadascii('path', fullfile(find_path_toolbox('rmnspec'),'data_pur'),'noprefetch');
%   peaks=monotone2peaks(dbpur.Stearicacid400.I,'mfilt',[1 10]);
%   [isbase,Ib] = nmrbaseline(dbpur.Stearicacid400.I,'mfilt',10,'prejection',peaks(2).preject);
%   figure,plot(dbpur.Stearicacid400.ppm,dbpur.Stearicacid400.I,'b-'), hold on,  plotequi(dbpur.Stearicacid400.ppm(isbase),dbpur.Stearicacid400.I(isbase,1),'r-'), plot(dbpur.Stearicacid400.ppm,peaks(2).If,'g')
%   legend({'spectrum' 'peaks' 'filtered'})
%
% Example: bad signal to noise ratio
%   if isempty(find_path_toolbox('rmnspec')), error('install first the toolbox rmnspec'), end
%   [dbpur,dbxpur] = nmrloadascii('path', fullfile(find_path_toolbox('rmnspec'),'data_pur'),'noprefetch');
%   peaks=monotone2peaks(dbpur.Stearicacid100.I,'mfilt',[1 27]);
%   [isbase,Ib] = nmrbaseline(dbpur.Stearicacid100.I,'mfilt',27,'prejection',peaks(2).preject);
%   figure,plot(dbpur.Stearicacid100.ppm,dbpur.Stearicacid100.I,'b-'), hold on,  plotequi(dbpur.Stearicacid100.ppm(~isbase),dbpur.Stearicacid100.I(~isbase,1),'r-'), plot(dbpur.Stearicacid100.ppm,peaks(2).If,'g')
%   legend({'spectrum' 'peaks' 'filtered'})
% 
% Example : bad S/N ratio
%   if isempty(find_path_toolbox('rmnspec')), error('install first the toolbox rmnspec'), end
%   [dbpur,dbxpur] = nmrloadascii('path', fullfile(find_path_toolbox('rmnspec'),'data_pur'),'noprefetch');
%   peaks=monotone2peaks(dbpur.MBOCA1250.I,'mfilt',[1 27]);
%   [isbase,Ib] = nmrbaseline(dbpur.MBOCA1250.I,'mfilt',27,'prejection',peaks(2).preject);
%   figure,plot(dbpur.MBOCA1250.ppm,dbpur.MBOCA1250.I,'b-'), hold on,  plotequi(dbpur.MBOCA1250.ppm(~isbase),dbpur.MBOCA1250.I(~isbase,1),'r-'), plot(dbpur.MBOCA1250.ppm,peaks(2).If,'g')
%   legend({'spectrum' 'peaks' 'filtered'})

% RMNSPEC v 0.1 - 11/09/12 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 06/10/12
% 
% History
% 29/09/12 updated algorithm to match monotone2peaks, updated example based on plotequi
% 05/10/12 add example to debug baseline problem
% 06/10/12 add examples (different cases of signal to noise ratio)
%
% Default
default = struct('mfilt',20,'prejection',5);

% argcheck
o = argcheck(varargin,default);
n = size(I,1);

% 
If = filtzero(I,o.mfilt);

% baseline procedure is based on the same algorithm used in monotone2peaks
[pp,lp,ap] = monotone(If,'+');
[pm,lm,am] = monotone(If,'-');
pp = pp+lp-1;
pall = intersect(pp,pm);% position of peaks
nall = length(pall);    % number of peaks
aall = NaN(nall,2);     % amplitude of peaks
[start,stop] = deal(NaN(nall,1)); % begining and end of peaks
[~,ip,iall] = intersect(pp,pall); aall(iall,1) = ap(ip);  start(iall) = pp(ip)-lp(ip)+1;
[~,im,iall] = intersect(pm,pall); aall(iall,2) = -am(im); stop(iall)  = pm(im)+lm(im)-1;
aall = max(aall,[],2);

% identify significant peaks
significant = find(aall>=prctile(aall,100-o.prejection));
isbase = true(n,1);
for i=1:length(significant)
    isbase(start(significant(i)):stop(significant(i))) = false;
end

% output
if nargout>0
    base = find(isbase);
    fit = polyfit(base,I(base),1);
    Ib = I - polyval(fit,(1:n)');
end

    