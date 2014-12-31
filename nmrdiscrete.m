function Ifuzzy = nmrdiscrete(ppm,I,varargin)
%NMRFUZZY transforms NRM spectra in 
%   Syntax : out = nmrfuzzy(ppm,I[,'property',value,...)
%              ppm: mx1 vector of chemical shift of NMR spectra
%                I: mxn vector of intensity of NRM sepctra
%   Property/value
%        'ppmfilt': filetring width (default = 0.01) (to be larger than ppmfilt)
%         'ppmtol': tolerance in ppm for peak reconstruction (default = 0.005)
%     'prejection': percentage of excluded peaks (default: 5=5%)
%           'ntst': alternative method to prejection (default = 30);
%
%    Output :
%          Ifuzzy : nxm vector of intensities corresponding to I with 0 almost everywhere execept at peaks
%
%    EXAMPLE:
%       if isempty(find_path_toolbox('rmnspec')), error('install first the toolbox rmnspec'), end
%       [dbpur,dbxpur] = nmrloadascii('path', fullfile(find_path_toolbox('rmnspec'),'data_pur'),'ppmstandards',[-0.05 0.15; 1.50 1.65; 7.2 7.32;],'ppmmin',-0.5,'ppmmax',12);
%       dbxpur.Id = nmrdiscrete(dbxpur.ppm,dbxpur.I)
%       figure, plot(dbxpur.ppm,dbxpur.Id)
%
% See also: nrlmoadascii, nmrcorr, nmrsubdb, nmrbaseline
%
% RMNSPEC v 0.1 - 14/09/12 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 17/09/12
% 
% History
% 16/09/12 release candidate (many fixes)
% 16/09/12 add vectorization column-wise, updated help, add example
% 17/09/12 add rejection optimizer

% Default
default = struct('ppmfilt',0.01,'ppmtol',0.005,'prejection',5,'ntst',1024);
%kappa = @(y1,y2) y2./(1+y1.^2).^1.5;


% argcheck
if nargin<2, error('2 arguments are required'), end
[m,n] = size(I);
if (size(ppm,1)~=m) || size(ppm,2)~=1, error('ppm must be a %dx1 array',m), end
o = argcheck(varargin,default);
Ifuzzy = zeros(m,n);
mfilt = max(1,ceil(o.ppmfilt./(ppm(2)-ppm(1))));
mtol  = max(1,ceil(o.ppmtol./(ppm(2)-ppm(1))));
mfilt = max(mtol,mfilt);
proba = linspace(0,100,o.ntst)';

% filtering
If = filtzero(I,mfilt);

% for each column of I
screen = ''; t0 = clock;
for k=1:n
    screen = dispb(screen,'nmrdiscrete %d/%d',k,n);
    % search monotone segments (peak = + followed by -)
    [pp,lp,ap] = monotone(If(:,k),'+');
    [pm,~,am] = monotone(If(:,k),'-');
    [pall,iall] = unique([pp+lp-1;pm]); % remove duplicate due to + and - process    
    aall = [ap;abs(am)];
    aall = aall(iall);
    amin1 = prctile(aall,100-o.prejection);    % user rejection
    criterion = ndf(proba,filtzero(prctile(aall,100-proba),10),2); % automatic rejection
    amin2 = prctile(aall,100-proba(5+find(criterion(5:end)<0,1,'first')));
    significant = find(aall>=min(amin1,amin2));
    significant = significant(I(pall(significant),k)>0);
    Iguess = zeros(m,1);
    Iguess(pall(significant)) = I(pall(significant),k);
    Ifuzzy(:,k) = conv(Iguess,ones(2*mtol+1,1),'same'); %>0;
    % check with: figure, plot(I(:,k)), hold on,  plot(Ifuzzy(:,k),'r-'), plot(If(:,k),'g-')
end
dispb(screen,'nmrdiscrete completed on %dx%d data in %0.4g s',m,n,etime(clock,t0));