function [gaussianpeak,outsum,outsumbaseline] = monotonepeakfit(p,varargin)
%MONOTONEPEAKFIT fits [x y] as a sum of Gaussian/Lorentzian peaks with positions previsously identified with monotonepeak
% Note that two methods of fit are always applied with prescribed centers (method 1) and movable centers (method 2)
% Results of both methods are returned.
% The simplified Gaussian kernel used in this function verifies: 
%                 area: 0.6006*sqrt(pi)*width
%   standard deviation: 0.6006/sqrt(2)*width
% The corresponding Lorenzian kernel shares the same value at half height
%
% SYNTAX: gaussianpeak = monotonepeakfit(p,'x',x,'y',y,'property1',value1...,'keyword')
%         [gaussianpeak,outsum,outsumwithbaseline] = monotonepeakfit(...)
%
% LIST OF PAIR PROPERTY/VALUE
%            p: npx1 structure array created with MONOTONEPEAK (with np = peak number)
%            x: vector of x values (default = [])
%            y: mx1 vector of y values
%  significant: percent to keep signifiant peaks (default = .95)  % search the first peaks >95% of weight
%minpointsinbaseline: (default 5)
%     shiftmax: maximum shift tolerated with method 2(default = [])
%      preject: rejection probability for baseline detection (default = 5);
%
% KEYWORDS
%   'baseline': fit and remove linear baseline 
%       'sort': to sort in an descend order of weigth attributed to gaussian peaks
% 'lorentzian': to fit lineshape with lorentzian peak
%  'endforced': to add last point of x to identify baseline
%  'keeporder': force Gaussians/Lorentzian to have increasing positions
%
% OUTPUTS (with the following notations: np = number of peaks, 2 = 2 fitting models)
%
%    gaussianpeak = np x 2 structure array with fields
%              rank: rank of gaussian peak (descendant order by weight)
%            weight: weight of gaussian peak 
%    relativeweight:relative weight of gaussian peak compaired to sum of weight
%             width: width of gaussian peak
%          position: position or center of gaussian peak
%          baseline: polyval function to define baseline of gaussian peak
%        isbaseline: index/point to define baseline of gaussian peak
%         xbaseline: point based on x to define baseline of gaussian peak
%            kernel: kernel function of gaussian peak
%                r2: coefficient of determination of gaussian peak
%             r2max: 2 maxima of coefficient of determination of gaussian peak
%             sigma: standard deviation of gaussian peak
%        expansion1: kernel function of all gaussian peaks without baseline (x, n = peak number, k = fitting models)
%        expansion2: kernel function of all gaussian peaks with baseline (x, n = peak number, k = fitting models)
%            window: is a structure with fields
%                   center: center of window 
%                    width: width of window (x(end)-x(1))
%                widthapha: width of window calculated with alpha
%                    alpha: probability of rejection in window (default=0.05)
%        nsignificantpeaks: number of significant peaks
%         significantproba: percentage to  keep signifiant peaks (variable 'significant', default = .95)
% 
%   outsum is an anonymous function (x,norder,imodel) giving the serial expansion
%               Sum i =1..norder weight(i) * kernel_i(x)      for imodel = 1 or 2
%   outsumwithbaseline is an anonymous function (x,norder,imodel) giving the serial expansion
%               (Sum i =1..norder weight(i) * kernel_i(x)) + baseline(x)       for imodel = 1 or 2
%
% See also monotonepeak, monotone
%
% BASIC example
%   g = @(x,position,width)  exp(-((x-position)./(0.6006.*width)) .^2)
%   x = linspace(0,1000,1000)';
%   y = 0.5*g(x,700,30) + 2*g(x,830,100);
%   p = monotonepeakfit; p.center = [690 850]; p.width = [10 20];
%   q = monotonepeakfit(p);
%   [q,~,model] = monotonepeakfit(p,'x',x,'y',y);
%   hp = plot(x,[y model(x,2,1) model(x,2,2)]); legend(hp,{'original' 'fit 1' 'fit with free center'},2)
%   
%
% ADVANCED EXAMPLE
% [dbpur,dbxpur] = nmrloadbase; 
% mol = {'Erucamide'}; % substance test Erucamide
% [~,isubs] = intersect(dbxpur.commonname,mol); % find index of Erucamide in databse
% ROI = dbpur.(mol{1}).gates; % extract RoI (nx4 array: ppmmin, ppmmax, buffer and weight)
% valid = ((dbxpur.ppm>=ROI(1,1))&(dbxpur.ppm<=ROI(1,2))); % take the 1st considered ROI
% x = dbxpur.ppm(valid);          
% y = dbxpur.I0(valid,isubs);                
% p = monotonepeak('x',x,'y',y,'array'); % find peaks
% [gaussianpeak,~,model] = monotonepeakfit(p,'x',x,'y',y,'baseline','sort'); % fitting
% subplot(211)
% hp = plot(x,[y model(x,size(gaussianpeak,1),1) model(x,size(gaussianpeak,1),2)]);
% hl = legend('measured','deconvoluted (s)','deconvoluted (p,s)');
% subplot(212)
% hp = plot(x,[y model(x,3,1) model(x,3,2)]);
%
%
% INRA\MS 2.1 - 24/03/2013 - Olivier Vitrac - rev. 23/11/2013
%
%
% TODO LIST
%   16/05/13 help on additional parameters, fix r2
%
% Revision history
% 25/03/13 add help and example
% 26/03/13 major update, remove cells, add outsum, add r2all, fix r2
% 27/03/13 add lorentzian peaks and relativeweight
%          fix help
% 29/03/13 add gaussianpeak.window for exact description of multiplet (center, width...)
% 09/04/13 add variable 'significant' to define rejection percent 
%          modify outsum = model without baseline
%                 outsumbaseline = model with baseline 
% 11/04/13 rename variables and add window.nsignificantpeaks,window.significantproba
% 16/04/13 add keyword 'endforced' to add last point of x to identify baseline
% 14/05/13 fix outputs, major help update
% 15/05/13 update help
% 23/11/13 accept mfilt=0

%% default
keyword = {'baseline','sort','lorentzian','endforced','keeporder'};
options = struct('display','iter','FunValCheck','on','MaxIter',1e3,'TolFun',1e-6,'TolX',1e-6);
default = struct('x',[],'y',[],'significant',.95,'options',options,'minpointsinbaseline',5,'shiftmax',[],'preject',5);
minumunrequiredfields = {'center','width','ibase'};
requiredfields = {'tail' 'wall' 'height' 'center' 'istart' 'start' 'stop' 'istop' 'width','ratioheight','ratiowidth','ibase'};

%% Kernels
% simplified Gaussian kernel
%   area: 0.6006*sqrt(pi)*width
%   standard deviation: 0.6006/sqrt(2)*width
% Lorenzian kernel shares the same value at half height
gaussiankernel = @(x,position,width)  exp(-((x-position)./(0.6006.*width)) .^2);
lorentziankernel =  @(x,position,width) 1./(1+((x-position)./(0.5.*width)).^2);

% argcheck
if nargin<1, gaussianpeak = cell2struct(repmat({[]},1,length(minumunrequiredfields)),minumunrequiredfields,2); return, end
if isempty(p), p = monotonepeak(varargin{:}); end
if ~isstruct(p) && all(cellfun(@(f) isfield(p,f),minumunrequiredfields)), error('p must be created with monotonepeak'), end
if ~all(cellfun(@(f) isfield(p,f),requiredfields)), dispf('WARNING:: object p has not been created with monotonepeak'), end
if length(p)==1 && length(p.center)>1, p = struct2structtab(p); end, n = length(p);
o = argcheck(varargin,default,keyword,'nostructexpand');
m = length(o.y);
if m == 0, error('y is empty'), end
if numel(o.y)~=m, error('y must be a vector'), end
if isempty(o.x), o.x = (1:m)'; end
if numel(o.x)~=m, error('x and y must be of the same size'); end
if isempty(o.shiftmax), o.shiftmax = min([p.width]); end
o.x = o.x(:); o.y = o.y(:);
sce = var(o.y)*(m-1);
penaltyscale = sqrt(sce)/n;
dx = median(diff(o.x))/4;
H = @(x) 1/2 * ( 1 + tanh(x/dx) ); % heaviside
% keyword 'lorentzian' if fitting with lorentzian peaks
if o.lorentzian, gaussiankernel = lorentziankernel; end

% remove baseline
if o.baseline
    isbaseline = true(m,1);
    for i=1:n
        isbaseline( (o.x>=p(i).start) & (o.x<=p(i).stop) ) = false;
    end
    if length(find(isbaseline))<o.minpointsinbaseline
        isbaseline([1 end])=true;
        isbaseline([p.ibase]) = true;
    end
    if o.endforced
        fluctuations = o.y;
        fluctuationsreject = prctile(fluctuations,o.preject);
        isbaseline(1:find(fluctuations>fluctuationsreject,1,'first')) = true;
        isbaseline(find(fluctuations>fluctuationsreject,1,'last'):m) = true;
    end
    b = polyfit(o.x(isbaseline),o.y(isbaseline),1);
    o.y = o.y - polyval(b,o.x);
else
    b = [0 0];
end

% Fitting
G = zeros(m,n);       % kernels (used by nested functions)
widthguess = cat(2,p.width)/4; %first guess of s
positionguess = cat(2,p.center); % first guess of p
[width,position,sigma]  = deal(NaN(2,n)); % 2 solutions - row-wise
[weight,cumweight,order,rank,relativeweight] = deal(NaN(n,2)); % 2 solutions - column-wise
[windowcenter,nsignificantpeaks] = deal(zeros(2,1));
critfit = zeros(2,1); % 2 solutions - column-wise
yfit    = NaN(m,2);   % 2 solutions - column-wise
% fit based on s only (peak positions set by [p.center])
warning off %#ok<WNOFF>
[width(1,:),critfit(1)] = fminsearch(@fitgaussianfixedpos,widthguess,o.options);
position(1,:) = positionguess;
% fit with positions and s free
[tmp,critfit(2)] = fminsearch(@fitgaussian,[positionguess;width(1,:)],o.options);
position(2,:) = tmp(1,:); width(2,:) = tmp(2,:);
% distance
warning on %#ok<WNON>
for k=1:2
    [yfit(:,k),weight(:,k)] = gaussian(o.x,position(k,:),width(k,:));
    [~,order(:,k)] = sort(weight(:,k),'descend');
    rank(order(:,k),k) = 1:n;
end

% control plot (to be removed)
% figure, plot(o.x,[o.y yfit]); hold on, plot(o.x(isbaseline),o.y(isbaseline),'r.')

% list of peaks
gaussianpeak = repmat(struct('rank',[],'weight',[],'relativeweight',[],...
                             'width',[],'sigma',[],'position',[],'baseline',[],'xbaseline',[],'isbaseline',[],...
                             'kernel',[],'r2',NaN,'r2max',[],...
                             'window',struct('center',[],'width',[],'widthalpha',[],'alpha',[],'nsignificantpeaks',[],'significantproba',[])),[n,2]);
alpha = 0.05; % probability of rejection
for k = 1:2
    relativeweight(:,k) = weight(:,k)/sum(weight(:,k),1);
    windowcenter(k) = position(k,:)*weight(:,k)/sum(weight(:,k),1);
    windowwidth = o.x(end)-o.x(1);
    windowwidthalpha = norminv([alpha/2 1-alpha/2],windowcenter(k),(0.6006*windowwidth)/sqrt(2));
    sigma(k,:) = (0.6006*width(k,:))/sqrt(2);
    for i=1:n
        gaussianpeak(i,k).kernel =  @(x) gaussiankernel(x,position(k,i),width(k,i));
        gaussianpeak(i,k).weight = weight(i,k);  
        gaussianpeak(i,k).relativeweight = relativeweight(i,k);
        gaussianpeak(i,k).rank = rank(i,k);
        gaussianpeak(i,k).r2max = 1 - critfit(k).^2/sce;
        if o.baseline
            gaussianpeak(i,k).baseline = @(x) polyval(b,x);
            gaussianpeak(i,k).isbaseline = isbaseline;
            gaussianpeak(i,k).xbaseline = o.x(isbaseline);
        end
        gaussianpeak(i,k).width = width(k,i);    
        gaussianpeak(i,k).position = position(k,i);        
        gaussianpeak(i,k).sigma = sigma(k,i);
        gaussianpeak(i,k).window.center = windowcenter(k);
        gaussianpeak(i,k).window.width = windowwidth;
        gaussianpeak(i,k).window.alpha = alpha;
        gaussianpeak(i,k).window.widthalpha = diff(windowwidthalpha);
    end
end

% sorting peaks if required
if o.sort || (nargout>1),
    % sort all models
    for k=1:2
        gaussianpeak(:,k) = gaussianpeak(order(:,k),k);
        position(k,:) = position(k,order(:,k));
        width(k,:) = width(k,order(:,k));
        weight(:,k) = weight(order(:,k),k);
        sigma(k,:) = sigma(k,order(:,k));
        cumweight(:,k) = cumsum(relativeweight(order(:,k),k),1);      
        nsignificantpeaks(k) = find(cumweight(:,k)>=o.significant,1,'first');
        windowcenter(k) = position(k,1:nsignificantpeaks(k))*weight(1:nsignificantpeaks(k),k)/sum(weight(1:nsignificantpeaks,k),1);       
    end
    % user override
    if o.keeporder
        [~,order] = sort([gaussianpeak(:,1).position]); gaussianpeak = gaussianpeak(order,:);
        position = position(:,order); width = width(:,order); weight = weight(order,:);        
    end
    % general expansion expression (for both models)
    expansion1 = @(x,n,k) gaussiankernel(repmat(x(:),1,n),position(k*ones(1,numel(x)),1:n),width(k*ones(1,numel(x)),1:n))*weight(1:n,k);
    expansion2 = @(x,n,k) gaussiankernel(repmat(x(:),1,n),position(k*ones(1,numel(x)),1:n),width(k*ones(1,numel(x)),1:n))*weight(1:n,k)+polyval(b,x(:));
    % use the expansion to refresh r2 and relativeweight with expansion order n
     for k=1:2
        for i=1:n
            gaussianpeak(i,k).r2 = 1 - norm(o.y-expansion2(o.x,i,k)).^2/sce;
            gaussianpeak(i,k).relativeweight =  cumweight(i,k);
            gaussianpeak(i,k).window.center = windowcenter(k);
            windowwidthalpha = norminv([alpha/2 1-alpha/2],windowcenter(k),(0.6006*windowwidth)/sqrt(2));
            gaussianpeak(i,k).window.widthalpha = diff(windowwidthalpha);
            gaussianpeak(i,k).window.significantproba = o.significant;
            gaussianpeak(i,k).window.nsignificantpeaks = nsignificantpeaks(k);
        end
     end
end

% additional outputs
if nargout>1, outsum = expansion1; end
if nargout>2, outsumbaseline = expansion2; end

%% Nested functions
    function err = fitgaussian(ps)
        err = norm(gaussian(o.x,ps(1,:),ps(2,:))-o.y) ...
              + norm((1-H(ps(2,:)))*penaltyscale) ...
              + norm((H( abs(ps(1,:)-[p.center])-o.shiftmax) )*penaltyscale);
    end

    function err = fitgaussianfixedpos(s)
        err = norm(gaussian(o.x,[p.center],s)-o.y) + norm((1-H(s))*penaltyscale);
    end

    function [y,weights] = gaussian(x,c,s)
        % base functions
        for j = 1:n
            G(:,j) = gaussiankernel(x,c(j),s(j));
        end
        % weights
        W = abs(G\o.y);
        % fit
        y = G*W;
        % additional out if required
        if nargout>1, weights = W; end
    end

end