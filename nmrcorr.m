function corr = nmrcorr(dbxref,dbxmix,varargin)
%nmrcorr calculates the correlation between two databases of NMR spectra (as loaded with NMRLOADASCII)
%     Syntax: out = nmrcorr(dbxref,dbxmix [,'property',value,keyword...])
%           dbxref: NMR spectra of reference substances (to be loaded with nmrloadascii)
%           dbxmix: NMR spectrum/spectra of mixture(s) (to be loaded with nmrloadascii)
%                   keep it empty when only pair comparisons are required
%     Property/value
%         'maxlag': max lag in ppm for inter/autocorrelation functions (default = 0.05)
%              'I': field used to calculate correlations (I or I0) (default = 'I')
%     Keyword
%            'all': force the correlation function to be retrieved
%            examples: corr = nmrcorr(dbxref,[],'all')
%                     corr = nmrcorr(dbxref,dbxmix,'all')
%
%     corr = structure with fields (by noting nref=dbxref.n, and npairs = nref(nref+1)/2
%            i: npairsx1 array of row subscripts corresponding to linear index ind
%            j: npairsx1 array of column subscripts corresponding to linear index ind
%          ind: npairsx1 array of linear index corresponding to (i,j)
%         ind2: npairsx1 array of linear index corresponding to (j,i)
%            C: npairsx1 array of maximum of pair correlations of dbxref (including autocorrelation ones)
%          rho: npairsx1 array of maximum of pair correlation coefficients of dbxref
%         rank: npairsx1 array of ranks of sorted pair correlation coefficients between of dbxref
%
%   If dbxmix is supplied, additional corr  fields include
%           Cm: nrefxdbxmix.n array of maximum of correlations between dbxref and dbxmix pairs
%         rhom: nrefxdbxmix.n array of maximum of correlation coefficients between dbxref and dbxmix pairs
%
%   If 'all' is used, additional corr fields include
%          lag: nlagsx1 
%          Cij: nlags x npairs array of pair correlation functions relative to dbxref
%         Cmij: nlags x nref x dbxmix.n array of pair correlation functions relative to dbxref and dbxmix
%
% See also: nrlmoadascii
%
%   Example (update it to your need, insufficiently tested to be considered for valid/safe for production)
% -------------------------------------
%     [dbpur,dbxpur] = nmrloadascii('odsfile','substance.ods','sheetname','substance','path','data_pur');
%     [dbmix,dbxmix] = nmrloadascii('odsfile','extract.ods','sheetname','extract','primarykey','reference','path','data_mixture');
%     corr = nmrcorr(dbxpur,dbxmix)
%     C  = zeros(dbxpur.n,dbxpur.n); C(corr.ind) = corr.C; C(corr.ind2) = corr.C; % build C with symmetry
%     clc
%     options = optimset('TolFun',1e-10,'TolX',1e-4,'MaxIter',1e6,'Display','iter','LargeScale','off');
%     minproba = 1e-2; % minimum probability to say that a compound is significant
%     minc2 = 0.3; % minimum probability to say that a compound is significant
%     for imix = 1:dbxmix.n
%         dispf('\n---------------------\nMIXTURE %d = %s\n---------------------',imix,dbxmix.material{imix})
%         Cm = corr.Cm(:,imix);
%         Constr = ones(1,dbxpur.n); %eye(dbxpur.n); Constr(1,:)=1;
%         constr = 1; %zeros(dbxpur.n,1); constr(1,1)=1;
%         p = lsqlin(C,Cm,[],[],Constr,constr,zeros(dbxpur.n,1),ones(dbxpur.n,1),zeros(dbxpur.n,1),options);
%         [psort,isort] = sort(p,'descend'); significant = 1:find(psort>minproba,1,'last');
%         if isempty(significant)
%             dispf('\nMIXTURE %d no result with linear solving (at p>%0.3g)',imix,minproba)
%         else
%             dispf('\nMIXTURE %d (%s)\n\t%d likelyest candidates based on linear solving',imix,dbxmix.material{imix},significant(end))
%             cellfun(@(n,c,p) dispf('%24s (%6s): p=%12.3e',n,c,p),dbxpur.commonname(isort(significant)),dbxpur.reference(isort(significant)),num2cell(psort(significant)))
%         end
%         maxc2 = 0.5;
%         [psort,isort] = sort(corr.rhom(:,imix),'descend'); significant = 1:find(psort>minc2,1,'last');
%         if isempty(significant)
%             dispf('\nMIXTURE %d no significant result with linear solving (at corr2>%0.3g)',imix,minc2)
%         else
%             dispf('\nMIXTURE %d (%s)\n\t%d likelyest candidates from pair correlations',imix,dbxmix.material{imix},significant(end))
%             cellfun(@(n,c,p) dispf('%24s (%6s): corr2=%12.3e',n,c,p),dbxpur.commonname(isort(significant)),dbxpur.reference(isort(significant)),num2cell(psort(significant)))
%         end
%     end

% RMNSPEC v 0.1 - 31/08/12 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev.11/09/12

% History
% 01/09/12 release candidate
% 02/09/12 full help, improved information display, add example
% 10/09/12 update example
% 11/09/12 improve help
% 11/09/12 implement nmrbaseline, add property 'I' to calculate correlations based on I or I0

%Default
default = struct('maxlag',0.05,'displaydelay',.5,'I','I');
keyword = 'all';

%argcheck
if ~nargin, error('1 argument is required'), end
if nargin<2, dbxmix = struct([]); end
o = argcheck(varargin,default,keyword,'case');
if ~isstruct(dbxref) || ~isfield(dbxref,'I') || ~isfield(dbxref,'ppm') || ~isfield(dbxref,'m')  || ~isfield(dbxref,'n')
    error('dbxref must be created with nmrloadascii')
end
if ~isempty(dbxmix) && (~isfield(dbxmix,'I') || ~isfield(dbxmix,'ppm') || ~isfield(dbxmix,'m')  || ~isfield(dbxmix,'n'))
    error('dbxmix must be created with nmrloadascii')
end

% OUPUTS preparation
sq2vec = @(Z) Z(Z>0);
maxlag = ceil(o.maxlag/dbxref.step);
ntot = dbxref.n*(dbxref.n+1)/2;
[i,j] = ndgrid(1:dbxref.n,1:dbxref.n);
corr = struct(...
    'i',sq2vec(tril(i,0)),...
    'j',sq2vec(tril(j,0)),...
    'ind',[],...
    'C',zeros(ntot,1),...
    'rho',zeros(ntot,1),...
    'rank',zeros(ntot,1) ...
        );
corr.ind = sub2ind(dbxref.n*[1 1],corr.i,corr.j);
corr.ind2 = sub2ind(dbxref.n*[1 1],corr.j,corr.i);
if o.all
    corr.lag = (-maxlag:maxlag)'*dbxref.step;
    corr.Cij = zeros(2*maxlag+1,ntot,'single');
end
    
% PAIR CORRELATIONS based on dbxref: C, Cij
V = var(dbxref.(o.I),1);
screen = ''; count = 0; t0 = clock; t1=t0; 
for j=1:dbxref.n
    for i=1:dbxref.n
        if i>=j % lower triangle plus main diagonal
            count = count + 1;
            if etime(clock,t1)>o.displaydelay
                screen = dispb(screen,'%dx%d pair correlations (lag<%d) %4.3g %% completed (remaining time %4.3g s)',...
                    dbxref.n,dbxref.n,maxlag+1,100*count/ntot,etime(clock,t0)*(ntot/count-1)); t1 = clock;
            end            
            c = xcorr(dbxref.(o.I)(:,i),dbxref.(o.I)(:,j),maxlag,'unbiased');
            corr.C(count) = max(c); % correlation
            corr.rho(count) = corr.C(count)/sqrt(V(i)*V(j)); % correlation coefficients
            if o.all
                corr.Cij(:,count) = single(c);
            end
        end
    end
end
dispb(screen,'NMRCORR: %dx%d pair correlations (lag<%d)  completed in %0.3g s',dbxref.n,dbxref.n,maxlag+1,etime(clock,t0));

% SORT rho values
[~,i] = sort(abs(corr.rho),'descend');
corr.rank(i) = 1:ntot;

% INTER CORRELATIONS between dbxref and dbxmix
if ~isempty(dbxmix)
    corr.Cm = zeros(dbxref.n,dbxmix.n);
    corr.rhom = zeros(dbxref.n,dbxmix.n);
    if o.all, corr.Cmij = zeros(2*maxlag+1,dbxref.n,dbxmix.n,'single'); end
    W = var(dbxmix.(o.I),1);
    ntot = dbxref.n * dbxmix.n;
    screen = ''; count = 0; t0 = clock; t1=t0; 
    for i=1:dbxref.n
        for j=1:dbxmix.n
            count = count + 1;
            if etime(clock,t1)>o.displaydelay
                screen = dispb(screen,'%dx%d inter correlations (lag<%d) %4.3g %% completed (remaining time %4.3g s)',...
                    dbxref.n,dbxmix.n,maxlag+1,100*count/ntot,etime(clock,t0)*(ntot/count-1)); t1 = clock;
            end
            c = xcorr(dbxref.(o.I)(:,i),dbxmix.(o.I)(:,j),maxlag,'unbiased');
            corr.Cm(count) = max(c); % correlation
            corr.rhom(count) = corr.Cm(count)/sqrt(V(i)*W(j)); % correlation coefficients
            if o.all, corr.Cmij(:,i,j) = single(c); end
        end
    end
    dispb(screen,'NMRCORR: %dx%d inter correlations (lag<%d) completed in %0.3g s',dbxref.n,dbxmix.n,maxlag+1,etime(clock,t0));
end


