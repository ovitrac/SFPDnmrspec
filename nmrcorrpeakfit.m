function corr = nmrcorrpeakfit(varargin)
%NMRCORRPEAKFIT calculate correlation coefficient between pattern (peak fitted) and NMR normalized spectrum (x,y)
% syntax: corr = nmrcorrpeakfit('x',x,'y',y,'xpattern',xfit,'ypattern',yfit,(property1',value1...,'sort')
%
% List of pair property/value%            
%            x: mx1 vector of x values (ppm unit)
%            y: mx1 vector of y values
%     xpattern: nx1 vector of x values of pattern ((ppm unit)(see dbfit of NMRLOADDBFIT)
%     ypattern: nx1 vector of y values of pattern (see dbfit of NMRLOADDBFIT)
%       nmatch: number of peak matched with pattern (default = 10)
%  
% KEYWORD
%       'sort': to sort in an descend order of correlation coefficient
%   'specific': to calculate correlation on spectifc zone (correspond with x of pattern), so that, it is not possible to use 'sort' keyword-> only one correlation coefficient will be calculated

% Outputs dbcorr: structure with size npx1 (np = number of matched areas )
%    corr.rank: npx1 rank of matching score (descendant order by correlation coeffcient)
%    corr.xfit: nx1 vector of x values of pattern
%    corr.yfit: nx1 vector of y values of pattern
%    corr.yref: nx1 vector of y values windowed  
%    corr.xref: nx1 vector of x values windowed
%corr.normcorr: array, correlation coefficient calculated between pattern and signal windowed in raw spectrum
%       corr.c: (2n-1)x4, correlation function [yrefyref yrefyfit yfityref yfit yfit]
%    corr.lags: (2n-1)x1, lags of correlation function
%    corr.cmax: array, maximum of correlation
%    corr.imax: array, index corresponding to maximum of correlation
%corr.czerolag: array, correlation at zerolag
%
% RMNSPEC v 0.1 - 22/04/2013 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev.
% history
% 30/04/13: add 'specific',[segments] as varargin so that user can define segment to be studied
% 06/05/13: add corr.cmax, corr.clagzero adn corr.imax
%% default
keyword = 'sort';
default = struct('x',[],'y',[],'xpattern',[],'ypattern',[],'nmatch',10,'specific',[]);
% Argcheck
o = argcheck(varargin,default,keyword);
m = length(o.y); n = length(o.ypattern); p = length(o.specific);
if m == 0, error('y is empty'), end
if numel(o.y)~=m, error('y must be a vector'), end
if numel(o.x)~=m, error('x and y must be the same size'); end
if n == 0, error('ypattern is empty'), end
if numel(o.ypattern)~=n, error('ypattern must be a vector'), end
if numel(o.xpattern)~=n, error('xpattern and ypattern must be the same size'); end

% main
step = o.x(2)-o.x(1);   
% y0 = o.y-mean(o.y); % y with zero mean
ypattern0 = o.ypattern-mean(o.ypattern); % ypattern with zero mean
maxlag = o.xpattern(end)-o.xpattern(1); % window of y pattern
if o.specific
    if p<2,
        disp('WARNING: user specific area is truncated, it was redefined with fitted signal')
        ROIvalid = ((o.x>=min(o.xpattern))&(o.x<=max(o.xpattern))); 
    else ROIvalid = ((o.x>=min(o.specific))&(o.x<=max(o.specific))); 
    end
    xrefspe = o.x(ROIvalid);
    yrefspe = o.y(ROIvalid); yrefspe0 = yrefspe-mean(yrefspe);  
    c = nmrcorrcoeff(yrefspe,o.ypattern);
    corr.xfit = o.xpattern;
    corr.yfit = c.y2;
    corr.xref = o.xpattern;
    corr.yref = c.y1;
    corr.normcorr = c.corrcoeff;
%     segment = (1:length(xrefspe)-length(o.ypattern))';  % segment to be added 
%     ytest = [ypattern0;segment*0]; % ypattern0 adjusted   
%     [c,lags] = xcorr([yrefspe0 ytest],'unbiased');
%     zerolag = length(xrefspe);
%     [cmax,imax] = max(c(:,2)); 
%     normcorr =  cmax./sqrt(c(zerolag,1)*c(zerolag,4));
%     corr.xfit = [o.xpattern;max(o.xpattern)+segment*step];
%     corr.yfit = [o.ypattern;segment*0];
%     corr.xref = xrefspe;
%     corr.yref = yrefspe;
%     corr.normcorr = normcorr; 
%     corr.c = c;
%     corr.lags = lags;
%     corr.cmax = cmax;
%     corr.imax = imax;
%     corr.czerolag = c(zerolag,2);
else
    [cglobal,lagsglobal] = xcorr(o.y,o.ypattern);  % corr global 
    valid = lagsglobal>=0; % truncate lags -> index in lagsglobal = index in whole spectrum
    p = monotonepeak('x',lagsglobal(valid),'y',cglobal(valid)/max(o.y),'mfilt',1/800*m,'sort','descend','array'); % search maxima in cglobal trapz(ppm,y)
    npeaktotest = min([o.nmatch size(p,1)]); %
    indpeak = [p(1:npeaktotest).icenter]; % index of maxima cglobal;
    xref = zeros(n,npeaktotest);
    yref = zeros(n,npeaktotest);
    yfit = zeros(n,npeaktotest);
    normcorr = zeros(npeaktotest,1);
%     cmax = zeros(npeaktotest,1);
%     imax = zeros(npeaktotest,1);
%     czerolag = zeros(npeaktotest,1);
%     c = cell(1,npeaktotest);
%     lags = zeros(2*n-1,npeaktotest);
    corr = repmat(struct('rank',[],'xfit',[],'yfit',[],'xref',[],'yref',[],'normcorr',[],'c',[],'lags',[],'cmax',[],'czerolag',[],'imax',[]),[npeaktotest,1]);
    for i = 1:npeaktotest
        xref(:,i) = (o.x(indpeak(i)):step:o.x(indpeak(i))+maxlag)'; % ppm to test (ppm(index)+window width)
        ytest = interp1(o.x,o.y,xref(:,i));
        c = nmrcorrcoeff(ytest,o.ypattern);
        yfit(:,i) = c.y2;
        yref(:,i) = c.y1;
        normcorr(i,1) = c.corrcoeff;
%         yref0 = yref(:,i)-mean(yref(:,i));
%         segment = (1:length(xref(:,i))-length(o.ypattern))';  % segment to be added  (1:(length(xref)-length(yfit))/2)';
%         ytest = [o.ypattern;segment*0]; % yfit adjusted  
%         ytest0 = ytest-mean(ytest);
%         [c{1,i},lags(:,i)] = xcorr([yref0 ytest0],'unbiased'); %,ceil(maxlag/step)
%         zerolag = length(xref(:,i)); czerolag(i,1) = c{1,i}(zerolag,2);
%         [cmax(i,1),imax(i,1)] = max(c{1,i}(:,2)); 
%         normcorr(i,1) =  cmax(i,1)./sqrt(c{1,i}(zerolag,1)*c{1,i}(zerolag,4));
    end
    [~,order] = sort(normcorr,'descend');
    rank(order) = 1:npeaktotest;

    % save in ouput structure
    for i = 1:npeaktotest    
        corr(i).rank = rank(i);
        corr(i).xfit = o.xpattern;
        corr(i).yfit = yfit(:,i);% ytest;
        corr(i).xref = xref(:,i);% corr(i).xref = xref(:,i);
        corr(i).yref = yref(:,i);
        corr(i).normcorr = normcorr(i,1); 
%         corr(i).c = c{1,i};
%         corr(i).lags = lags(:,i);
%         corr(i).cmax = cmax(i,1);
%         corr(i).imax = imax(i,1);
%         corr(i).czerolag = czerolag(i,1);
    end
    if o.sort
        corr(:) = corr(order(:));
    end
end


%% plots if no output
if ~nargout
    if o.specific
        figname = 'correlation';
        hfig = figure; formatfig(hfig,'figname',figname,'paperposition',[0.3387    0.7975   29.0000   19.3890],'paperorientation','landscape');
%         hs = subplots([1 1],1,0.05,0);
%         subplot(hs(1)), plot([corr.yref/var(corr.yref) corr.yfit/var(corr.yfit)])
        plot([corr.yref/var(corr.yref) corr.yfit/var(corr.yfit)])
        set(gca,'xdir','reverse')
        title(['corr. coeff.: ' num2str(corr.normcorr)])
        formatax(gca,'fontsize',8,'xlim',[0 length(corr.xref)],... hs(1)
                 'xtick',linspace(0,length(corr.xref),3),...
                 'xticklabel',formatsci(linspace(min(corr.xref),max(corr.xref),3),'economy'))
        hl = legend({'raw signal' 'fitted signal'}); set(hl,'fontsize',5,'box','off','Color','none')
%         subplot(hs(2)), plot(corr.lags,corr.c), xlim([-length(corr.xref)+1 length(corr.xref)-1])  
%         formatax(hs(2),'fontsize',8)
    else
        figname = 'correlation';
        hfig = figure; formatfig(hfig,'figname',figname,'paperposition',[0.3387    0.7975   29.0000   19.3890],'paperorientation','landscape');
        hs = subplots(1,[0.05 0.95],0,0.05);
        nrows = ceil(sqrt(npeaktotest));
        ncols = ceil(npeaktotest/nrows);
        signalinfo = sprintf('Fitted signal\nPosition: [%0.3g;%0.3g] ppm',min(corr(1,1).xfit),max(corr(1,1).xfit));
        subplot(hs(1)), text(.5,0,signalinfo,'fontsize',10,'verticalalignment','bottom','horizontalalignment','center')
        subplot(hs(2)), hs1 = subplots(ones(1,ncols),ones(1,nrows),0.04,0.06);
        for i = 1:npeaktotest
%             subplot(hs1(i)), hs2 = subplots([1 1],1,0.02,0);
%             subplot(hs2(1)), plot([corr(i,1).yref/var(corr(i,1).yref) corr(i,1).yfit/var(corr(i,1).yfit)])
            subplot(hs1(i)), plot([corr(i,1).yref/var(corr(i,1).yref) corr(i,1).yfit/var(corr(i,1).yfit)])
            set(gca,'xdir','reverse')
            title(['corr. coeff.: ' num2str(corr(i,1).normcorr)])
            formatax(gca,'fontsize',8,'xlim',[0 n],... hs2(1)
                     'xtick',linspace(0,n,3),...
                     'xticklabel',formatsci(linspace(min(corr(i,1).xref),max(corr(i,1).xref),3),'economy'))
            hl = legend({'raw signal' 'fitted signal'}); set(hl,'fontsize',5,'box','off','Color','none')
%             subplot(hs2(2)), plot(corr(i,1).lags,corr(i,1).c(:,2)), xlim([-n+1 n-1])  
%             formatax(hs2(2),'fontsize',8)
        end
        delete(hs1(i+1:end))
        set(hs(1),'visible','off')
    end
end