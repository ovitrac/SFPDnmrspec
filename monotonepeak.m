function peaks = monotonepeak(varargin)
%MONOTONEPEAK find peak position based on monotone and filtzero
%   SYNTAX: p=monotonepeak('property1',value1,'property2',value2,...'keyword')
%
%   LIST OF PAIR PROPERTY/VALUE
%        x: vector of x values (default = [])
%        y: mx1 vector of y values (mandatory)
%    mfilt: filter bandwidth, see filtzero (default = m/100)
%          several values can be applied if peaks must be found with different signal-to-noise ratios
%mfiltsort: sorting methods of mfilt values (default = 'descend')
%     sort: 'ascend' or 'descend' (default ='' no sorting)
%     xtol: tolerance on x to identifiy similar peaks when several mfilt values are applied (default = [])
%     itol: index tolerance on x (itol is overridden by xtol if defined) to identify similar peaks (default = 10);
%  maxpeak: Inf maximum of peaks to return ('array' is forced)
%
%   KEYWORDS
%           'array' is equivalent to p = struct2structtab(p)
%       'keeporder' force peaks to have increasing centers
%
%   OUTPUTS
%   p = structure with fields (assume that the number of peaks is np)
%           level: 1 except if several mfilt values are used. In this case, peak found with mfilt(level), sorted according to mfiltsort
%            tail: 1xnp cell array of 'left' or 'right' for left or right tailed respectively
%            wall: 1xnp cell array of 'left' or 'right' 
%          height: peak height (according to monotone segment, 1xnp array)
%       absheight: peak height (corresponding to y values, 1xnp array)
%          center: peak position (1xnp array)
%           start: start position (1xnp array)
%            stop: stop position (1xnp array)
%           width: width (1xnp array)
%     ratioheight: vertical shape factor >=1 (1xnp array)
%      ratiowidth: horizontal shape factor >=1 (1xnp array)
%   If 'array' is used p is a npx1 structure arra
%
%
%   See also: monotone, filtzero, monotonepeak, monotone2peaks (toolbox: RMNspec)

% INRA\MS 2.1 - 20/03/2013 - Olivier Vitrac - rev. 23/05/2013

% Revision history
% 21/03/2013 add istart, istop
% 31/03/2013 fix double for x and y, add graphical output if required
% 11/04/2013 implement multiple mfilt values
% 07/05/2013 add absheight
% 14/05/2013 default value of 'sortby' set to absheight
% 14/05/2013 fix ibestbase (when only one point available)
% 14/05/2013 add 'keeporder', 'maxpeak'
% 23/05/2013 fix empty peaks when no peak is found
%            add 'zero' for monotone as varargin

% default
keyword = {'array','keeporder'};
default = struct('x',[],'y',[],'mfilt',[],'sort','','sortby','absheight','itol',10,'xtol',[],'mfiltsort','descend','maxpeak',Inf,'zero',[]);

mfiltratio_default = 1/100;

% arg check
o = argcheck(varargin,default,keyword);
o.y = o.y(:); o.x = o.x(:);
m = length(o.y);
if m == 0, error('y is empty'), end
if numel(o.y)~=m, error('y must be a vector'), end
if isempty(o.x), o.x = (1:m)'; end
if numel(o.x)~=m, error('x and y must be of the same size'); end
if isempty(o.mfilt), o.mfilt = mfiltratio_default*m; end
if ~isa(o.x,'double'), o.x = double(o.x); end
if ~isa(o.y,'double'), o.y = double(o.y); end
if ~nargout, o.array = true; o.sort = 'descend'; end
if any(o.xtol), o.itol = max(round(o.xtol/mean(diff(o.x))),1); end
o.mfilt = max(0,round(o.mfilt(:)));
nfilt = length(o.mfilt);

if nfilt==1 % one single mfilt value
    
    % filtering
    yf = filtzero(o.y,o.mfilt);
    
    % monotone + and -
    
    [pp,lp,ap] = monotone(yf,'+',o.zero);
    [pm,lm,am] = monotone(yf,'-',o.zero); am = -am;
    [p,ip,im] = intersect(pp+lp-1,pm);
    lpx = o.x(pp+lp-1)-o.x(pp);
    lmx = o.x(pm+lm-1)-o.x(pm);
    tail = {'left' 'right'};
    peaks = struct(...
        'level', ones(1,length(p)),... by default such that 1 means associated to mfilt(1)
        'tail',{tail((lpx(ip)<lmx(im))+1)},...
        'wall',{tail((ap(ip)<am(im))+1)},...
        'height', max(ap(ip),am(im))',...
        'absheight',o.y(p'),...
        'icenter',p',...
        'center',o.x(p)',...
        'istart', pp(ip)',...
        'start',  o.x(pp(ip))',...
        'istop', pm(im)'+lm(im)'-1,...
        'stop',   o.x(pm(im)+lm(im)-1)',...
        'width',  lpx(ip)'+lmx(im)',...
        'ratioheight',  max(ap(ip),am(im))'./min(ap(ip),am(im))',...
        'ratiowidth',  max(lpx(ip),lmx(im))'./min(lpx(ip),lmx(im))',...
        'ibase',[] ...
        );
    if ~isempty(peaks.istart)
        base = [peaks.istart;peaks.istop];
        [~,ibestbase] = min(o.y(base),[],1);
        ibestbase = ibestbase(1:size(base,2)); % added 14/05/13
        peaks.ibase = base(sub2ind(size(base),ibestbase,(1:size(base,2))));
    end
    % sort if requested
    if ~isempty(o.sort)
        [~,ind] = sort(peaks.(o.sortby),o.sort);
        for f = fieldnames(peaks)'
            peaks.(f{1}) = peaks.(f{1})(ind);
        end
    end
    % array forced
    if (o.array || o.keeporder || ~isinf(o.maxpeak))
        if ~isempty(peaks.istart)
            peaks = struct2structtab(peaks);
        else
            f = fieldnames(peaks);
            peaks = repmat(cell2struct(repmat({[]},length(f),1),f),0,1);
        end
    end
    
else %% if several mfilt values
    o.mfilt = sort(o.mfilt,o.mfiltsort);
    peaks = monotonepeak('array',true,'mfilt',o.mfilt(1),o);
    for i=2:nfilt %for each additional mfilt value
        ptmp = monotonepeak('array',true,'mfilt',o.mfilt(i),o);
        for j=1:length(ptmp) % for each new peak, we check that no similar peak exists in pref
            if min(abs(ptmp(j).icenter-[peaks.icenter]))>o.itol % the center of the new peak is sufficiently different from previous ones
                peaks(end+1) = ptmp(j); %#ok<AGROW>
                peaks(end).level = i;
            end
        end % next new peak
    end % next mfilt value
    if ~isempty(o.sort)
        [~,ind] = sort([peaks.(o.sortby)],o.sort);
        peaks = peaks(ind);
    end    
end

%% user override - added 14/05/13
% apply maxpeak
if  ~isinf(o.maxpeak), peaks = peaks(1:min(o.maxpeak,length(peaks))); end
% apply keeporder
if o.keeporder, [~,order] = sort([peaks.center]); peaks = peaks(order); end

%% plots if no output
if ~nargout
    hold on
    n = length(peaks);
    col = flipud(cbrewer('qual','Set2',n));
    plot(o.x,o.y,'k','linestyle','-','linewidth',.5)
    for i=1:n
        ind = peaks(i).istart:peaks(i).istop;
        plot(o.x(ind),o.y(ind),'linestyle','-','linewidth',2,'color',col(i,:))
        text(peaks(i).center,o.y(peaks(i).icenter),sprintf('%d',i),'VerticalAlignment','bottom','HorizontalALignment','center','fontsize',8,'fontweight','bold');
        drawnow, if i ==1, y0 = min(ylim); end
        line([1;1]*peaks(i).center,y0+[0;peaks(i).(o.sortby)],'linestyle','-','color',col(i,:),'linewidth',2)
        plot(peaks(i).center,y0+peaks(i).(o.sortby),'marker','o','markersize',6,'markerfacecolor',col(i,:),'markeredgecolor',col(i,:),'linestyle','none')
        line([peaks(i).start;peaks(i).center;peaks(i).stop],y0+[0;peaks(i).(o.sortby);0],'linestyle','-','color',col(i,:),'linewidth',1)
    end
end