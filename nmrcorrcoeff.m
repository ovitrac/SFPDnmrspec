function [corrcoeff,sumyrefyout,corr] = nmrcorrcoeff(yref,y,lags,varargin)
%NMRCORRCOEFF calcualtes correlation coefficient 
% syntax: [corrcoeff,sumy,corr] = nmrcorrcoeff(y1,y2,[])
%         [corrcoeff,sumy,corr] = nmrcorrcoeff(y1,y2,lags,...)
% INPUTS
%           yref: m1x1 vector of y values (yref is movable relative to y)
%              y: m2x1 vector of y values (y is fixed)
%           lags: array, lags to calculate correlation (in point unit) 
%      gateproba: probability to set gate to calculate corr.coeff. (default = .001)
%         method: extrapolation method (to set gate) (default = cubic)
%          scale: method to scale yref and y: 'energy', 'linear', 'variance' (default = 'energy')
% KEYWORD:
%       unbiased: 
%           sort: sort by raw correlation coefficients in a descend order
% OUPUTS 
%      corrcoeff: 1x3 array, correlation coefficients calculated by 3 methods 
%                 [raw all window, gate with n = lenght of gate, gate inf]
%       sumyrefy: array, value corresponding to sum(yref*y,1)
%           CORR: structure with fields
%           yref: m+npadx1 vector, windowed yref
%       yrefnorm: m+npadx1 vector, windowed and normalized yref for 1 between 3 modes ('energy', 'linear' or 'variance'), gate not applied
%              y: m+npadx1 vector, windowed y
%          ynorm: m+npadx1 vector, windowed and normalized y for 1 between 3 modes ('energy', 'linear' or 'variance'), gate not applied
%            rho: equation to calculate corr.coeff
%        typepad: type od padding according to lags
%         corr.m: array, length of yref and y after zeropadding 
%      corr.npad: array, padding to be added according to lags
%      corr.lags: array, corresponding lags  
%  corr.xgatemin: array, right border of gate
%  corr.xgatemax: array, left border of gate
%
% CASE 1: yref and y are vectors (yref (m1xn1) and y (m2xn2) with n1=1 & n2=1, lags (mlagsx1) with mlags>=1)
% INPUTS
%           yref: m1x1 vector of y values (yref is movable relative to y)
%              y: m2x1 vector of y values (y is fixed)
%           lags: mlagsx1 vector, lags to calculate correlation (in point unit)  
% OUTPUTS 
%      corrcoeff: mlagsx3 vector, correlation coefficients for all lags calculated by 3 methods
%       sumyrefy: mlagsx1 vector, sum(yref*y,1) for all lags
%           CORR: structure with size mlagsx1
%       
% CASE 2: yref and y are matrix (yref (m1xn1) and y (m2xn2) with n1>1 or n2>1, lags (1xnlags) with nlags = n1xn2)
% INPUTS
%           yref: m1xn1 vector of y values (yref is movable relative to y)
%              y: m2xn2 vector of y values (y is fixed)
%           lags: 1xnlags vector, lags to calculate correlation (in point unit)  
% OUTPUTS 
%      corrcoeff: n1xn2x3 matrix, correlation coefficients for all pairs of yref and y calculated by 3 methods
%       sumyrefy: n1xn2 vector, sum(yref*y,1) for all pairs of yref and y
%           CORR: structure with size n1xn2
%
% EXAMPLE
% dbfit = nmrloaddbfit('path',fullfile(find_path_toolbox('rmnspec'),'data_pur'),'dbname','dbfit.mat');
% [~,dbxpur,~,~,~,~,dbxmix] = nmrloadbase;
% fitmol = fieldnames(dbfit); 
% ppm =dbxpur.ppm;    
% imix=98; %SFPDPP3
% y = dbxmix.I(:,imix);
% ifit = 15; iROI = 2; % Erucamide
% xfit = dbfit.(fitmol{ifit}).all(iROI,1).ppm;
% yfit = dbfit.(fitmol{ifit}).all(iROI,1).Ifit;
% valid = ((ppm>=min(xfit))&(ppm<=max(xfit)));
% yreal = y(valid);
% [corrcoeff,sumy,corr] = nmrcorrcoeff(yfit,yreal,[],'unbiased');
% [corrcoeff,sumy,corr] = nmrcorrcoeff(yfit,yreal,62);

% see also: nmrfindlags
% RMNSPEC v 0.1 - 06/05/2013 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev.23/05/2013
% history
% 21/05/13 change output structure ([corrcoeff,sumyrefyout,corr.(fields)) instead of [corrcoeff,sumy1y2] 
% 23/05/13 add gate to calculate correlation coefficient
%          update help
% 24/05/13 add scale for intensity normalization by (inear regression, energy and variance)
% 27/05/13 fix equation of rho2
% 26/11/14 cubic method replanced by 'pchip' (remove Warning: INTERP1(...,'CUBIC') will change in a future release. Use INTERP1(...,'PCHIP') instead. )

% default
keyword = {'unbiased','sort'};
options = struct('display','iter','FunValCheck','on','MaxIter',1e3,'TolFun',1e-6,'TolX',1e-6);
default = struct('ploton',~nargout,'axes',gca,'gateproba',.001,'method','pchip','options',options,'scale','energy');

% definition of correlation coefficient r2
rho2 = @(y1,y2,n) max(((sum(y1.*y2,1)-sum(y1,1)*sum(y2,1)/n)./sqrt((sum(y1.^2,1)-sum(y1,1).^2/n)*(sum(y2.^2,1)-sum(y2,1).^2/n))),0); 

% argcheck
o = argcheck(varargin,default,keyword);
if nargin<2, error('2 arguments are required'), end
if nargin<3, lags = []; end
if isvector(yref), yref = yref(:); end
if isvector(y), y = y(:); end
[m1,n1] = size(yref); [m2,n2] = size(y); m = max(m1,m2);
[mlags,nlags] = size(lags); 

%--------------------------------- RECURSION ---------------------------------------------------------------------
% CASE 1: y1 and y2 are vectors (y1 (m1x1) and y2 (m2x1)) and lags (mlagsx1) with mlags>1
if n1==1 && n2==1 && mlags>1
    % definition of subplot
    if o.ploton, nrows = ceil(sqrt(mlags)); ncols = ceil(mlags/nrows);
        hs = subplots(ones(1,ncols),ones(1,nrows),0.04,0.08);
    end
    corrcoeff = zeros(mlags,3);
    sumyrefy = zeros(mlags,1);
    for i = 1:mlags
        if o.ploton, 
            nmrcorrcoeff(yref,y,lags(i),'axes',subplot(hs(i)),'gateproba',o.gateproba,'scale',o.scale)
        else  tmp = [];
            [corrcoeff(i,:),sumyrefy(i,1),tmp] = nmrcorrcoeff(yref,y,lags(i),'gateproba',o.gateproba,'scale',o.scale);
            if i ==1, corr = repmat(tmp,[mlags,1]);
            else corr(i) = tmp; end
        end
    end
    if nargout>1, sumyrefyout = sumyrefy; end
    if nargout>2, corr = corr; end
    if o.sort && ~o.ploton, [~,order] = sort(corrcoeff(:,1),'descend'); % sort by descend order of raw corr. coeff.
        corrcoeff = corrcoeff(order,:);
        sumyrefyout = sumyrefyout(order);
        corr(:) = corr(order); 
    end
    return
% CASE 2: y1 and y2 are matrix (y1 (m1xn1) and y2 (m2xn2) with n1>1 or n2>1, lags (1xnlags) with nlags = n1xn2)
elseif (n1>1 || n2>1) && mlags==1 
    if nlags~=n1*n2, error('please check lags supplied'), end
    % definition of suplot
    if o.ploton, nrows = ceil(sqrt(nlags)); ncols = ceil(nlags/nrows);
        hs = subplots(ones(1,ncols),ones(1,nrows),0.04,0.08);
    end
    l = 0;
    corrcoeff = zeros(n1,n2,3);
    sumyrefy = zeros(n1,n2);
    for i = 1:n1
        for j = 1:n2
            l=l+1;
            if o.ploton, 
                nmrcorrcoeff(yref(:,i),y(:,j),lags(l),'axes',subplot(hs(l)),'gateproba',o.gateproba,'scale',o.scale)
            else  tmp = [];
                [tmpcorr,sumyrefy(i,j),tmp] = nmrcorrcoeff(yref(:,i),y(:,j),lags(l),'gateproba',o.gateproba,'scale',o.scale);
                corrcoeff(i,j,:) = permute(tmpcorr,[1 3 2]);
                if (i==1) && (j==1), corr = repmat(tmp,n1,n2);
                else corr(i,j) = tmp; end
            end
        end
    end    
    if nargout>1, sumyrefyout = sumyrefy; end
    if nargout>2, corr = corr; end
    return
elseif n1>1 || n2>1 || mlags>1
    error('not implemented case')
end    

%---------------------------------------------------------------------------------------------
% STEP 1: fix size by right padding with zero
yref = [yref;zeros(m-m1,1)];
y = [y;zeros(m-m2,1)];  

% STEP 2: find lags if not supplied
if isempty(lags), 
    lags = nmrfindlags(yref,y,'nmax',1); % lags corresponding to maxima of correlation
end

% STEP 3: gate definition
% the gate is a weighted function with 1 in the region of interest of the yref and 0 outside
x = (1:m)';
Iref = cumtrapz(x,yref); Iref = Iref/Iref(end);
[~,imin] = min(abs(Iref-o.gateproba));
[~,imax] = min(abs(Iref-(1-o.gateproba)));
xmin = round(fminsearch(@(x0) abs(o.gateproba-interp1(x,Iref,x0,o.method,1)),x(imin),o.options));
xmax = round(fminsearch(@(x0) abs((1-o.gateproba)-interp1(x,Iref,x0,o.method,0)),x(imax),o.options));
gate = zeros(m,1); gate((x>=xmin)&(x<=xmax)) = 1;

% STEP 4: find pad to be added
npad = abs(lags);
pad = zeros(npad,1); 
if lags>0
    type = 'padrefright';
    gate = [gate;pad]; yref = [yref;pad]; y = [pad;y];
else
    type = 'padrefleft';
    gate = [pad;gate]; yref = [pad;yref]; y = [y;pad]; 
end

% STEP 5: correlation coefficient calculation
yrefgate = yref.*gate; ygate = y.*gate;
corrcoeff = zeros(1,3);
if max(yref) <= 0, corrcoeff(1,1) = rho2(yref,y/max(y),m+npad); 
elseif max(y) <= 0, corrcoeff(1,1) = rho2(yref/max(yref),y,m+npad); 
else corrcoeff(1,1) = rho2(yref/max(yref),y/max(y),m+npad);
end % raw corr.coeff (all window, no gate)
if max(yrefgate)<= 0 && max(ygate) <= 0
    corrcoeff(1,2) = rho2(yrefgate(yrefgate~=0),ygate(yrefgate~=0),xmax-xmin);  % corr.coeff with gate lenght 
    corrcoeff(1,3) = rho2(yrefgate,ygate,m+npad); % corr.coeff with gate and window lenght
elseif max(yrefgate) <= 0
    corrcoeff(1,2) = rho2(yrefgate(yrefgate~=0),ygate(yrefgate~=0)/max(ygate),xmax-xmin);
    corrcoeff(1,3) = rho2(yrefgate,ygate/max(ygate),m+npad);
elseif max(ygate) <= 0
    corrcoeff(1,2) = rho2(yrefgate(yrefgate~=0)/max(yrefgate),ygate(yrefgate~=0),xmax-xmin); % corr.coeff with gate lenght 
    corrcoeff(1,3) = rho2(yrefgate/max(yrefgate),ygate,m+npad); % corr.coeff with gate and window lenght
else
    corrcoeff(1,2) = rho2(yrefgate(yrefgate~=0)/max(yrefgate),ygate(yrefgate~=0)/max(ygate),xmax-xmin); % corr.coeff with gate lenght 
    corrcoeff(1,3) = rho2(yrefgate/max(yrefgate),ygate/max(ygate),m+npad); % corr.coeff with gate and window lenght 
end

if o.unbiased, sumyrefy = sum(yref.*y,1)/(m+npad);
else sumyrefy = sum(yref.*y,1)/m;
end

% save in ouput structure
switch lower(o.scale)
    case 'energy'
        yrefnorm = yref/sum(yref.^2,1); ynorm = y/sum(y.^2,1);
    case 'linear'
        P = polyfit(yref,y,1);
        yrefnorm = yref; ynorm = (y - P(2))/P(1); 
    case 'variance'
        yrefnorm = yref/var(yref); ynorm = y/var(y);
end
if nargout>1, sumyrefyout = sumyrefy; end
if nargout>2
    corr = struct('yref',[],'yrefnorm',[],'y',[],'ynorm',[],'lags',[],'xgatemin',[],'xgatemax',[],'rho2',[],'npad',[],'typepad',[],'m',[]);
    corr.yref = yref;
    corr.y = y;
    corr.lags = lags;
    corr.rho2 = rho2;
    corr.npad = npad;
    corr.typepad = type;
    corr.m = m;
    corr.xgatemin = xmin; corr.xgatemax = xmax;
    corr.yrefnorm = yrefnorm;
    corr.ynorm = ynorm;
end

%% plots if no output
if o.ploton     
    plot([yrefnorm ynorm]), legend({'yref', 'yreal'})
    title(sprintf('\\rho^2_{raw}=%0.4g; \\rho^2_{gate}=%0.4g; \\rho^2_{gateinf}=%0.4g',corrcoeff))  
end