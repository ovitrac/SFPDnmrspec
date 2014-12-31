function plotnmrlsqnonneg(dataplot,dbout,varargin)
%PLOTNMRLSQNONNEG plots result of deconvolution of NMR spectra calculated by NMRLSQNONNEG (C or Crelative as outputs)
% SYNTAXE  figure, plotnmrlasnonneg(W,ind)
% INPTUS
%    dataplot: structure of data (4th output of NMRLSQNONNEG)
%       dbout: structure of data (3th output of NMRLSQNONNEG)
% VARARGIN
%            'class': n x 1 array, rang of susbtances as used in nmrlsqnonneg
%          'figname': name of figure for save (default='deconvolution)
%           'parent': handle to plot
%    'paperposition': position of figure, (default = [])
% 'paperorientation': orientation of paper (default='landscape')
%        'papertype': type of paper (A3,A4...), (default = 'A4')
%           'colmol': color scale for susbtances (default=[])
%           'marker': n x 1 cell for marker for substances (default=suite de lettre alphabetic [a,b...])s
%       'markersize': numeric, size of marker (default=8)
%        'fontsize1': fontsize for axes (default=12)
%        'fontsize2': fontsize for legends (default=10)
%           'ylabel': text for label of y axis (default='concentration')
%         'nmolplot': numeric, number of substances to be plotted (default = 5)
%         'nplotfit': numeric, number of susbtances to be fitted and plotted (default=3)
%           'yscale': 'linear' or 'log' for y axis scale
%             'ylim': lower limite of y axis (default=0) (if log scale, ylim>0)
%          'rhosort': n x 1 array of correlation coefficient between susbtances and mixture spectrum
%   'margincolorbar': numeric (<1) for set axis subplot to plot colorbar (default = 0.1)
%      'margingraph': numeric (<1) for set axis subplot to plot tree graph of scenarios (default = empty)
% KEYWORD 
%          'plotfit': on to plot fitting data
%
% See also: nmrlsqnonneg, createcombination
%
% RMNSPEC v 0.1 - 05/12/2013 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 23/09/14
% History revision
% 23/09/2014 : main modification : fig alike figs 11 and 12 in NMR article submited in I&EC
%
% default
default = struct('class',[],'parent',[],'figname','deconvolution','paperposition',[],'paperorientation','landscape','papertype','A4','colmol',[],'marker','','markersize',8,'fontsize1',12,...
                 'fontsize2',10,'ylabel','concentration','nmolplot',5,'nplotfit',3,'yscale','log','ylim',1e-3,'rhosort',[],'margincolorbar',0.1,'margingraph',[],'rhothreshold',0.75);
keyword = 'plotfit';
% argcheck
o = argcheck(varargin,default,keyword);
if nargin < 2, error('2 arguments are required'), end
if ~isstruct(dataplot) || ~isfield(dataplot,'xplot') || ~isfield(dataplot,'ninsert') || ~isfield(dataplot,'Wtabout') || ~isfield(dataplot,'Wreltabout') || ~isfield(dataplot,'Ctabout') || ~isfield(dataplot,'Caveragetab') || ~isfield(dataplot,'presencematrix') || ~isfield(dataplot,'graph')    
    error('DATAPLOT must be created by NMRLSQNONNEG (4th output element)')
end
if ~isstruct(dbout) || ~isfield(dbout,'combinations') || ~isfield(dbout,'valuesmin') || ~isfield(dbout,'valuesmax') || ~isfield(dbout,'indmol') || ~isfield(dbout,'T') || ~isfield(dbout,'N')    
    error('DBOUT must be created by NMRLSQNONNEG (3th output element)')
end
if ~isempty(o.class) && size(o.class,1) ~= size(dataplot.Caveragetab,1), error('Size of class array and dataplot.Caveragetab must be equal (in row)'), end
if isempty(o.class), o.class = 1:size(dataplot.Caveragetab,1); end

% plot
% general axes for plot
if isempty(o.parent) 
    if ~isempty(o.figname), figname = o.figname; else figname = 'deconvolution'; end
    hfig = figure; formatfig(hfig,'figname',figname,'paperposition',o.paperposition,'paperorientation',o.paperorientation,'papertype',o.papertype)
end
hs = subplots([.48 .52],1,.01,0,'strict');
if ~isempty(o.colmol), colmol= o.colmol; else colmol = cbrewer('qual','Dark2',size(dataplot.Caveragetab,1)); end
if ~isempty(o.marker), marker = o.marker; else marker = cellfun(@(m) sprintf('%s',char(96+m)),num2cell(1:size(dataplot.Caveragetab,1)+1),'UniformOutput',false); end
% tree size according to bumber of branches
if ~isempty(o.margingraph), margingraph = o.margingraph; 
else
    nbranch = length(dbout.pathsN);
    margingraph = 1/(nbranch+2);
end
% horizontal bar according to rhothreshold
ifalsepos = length(o.rhosort(o.rhosort > o.rhothreshold));
subplot(hs(1)), hold on
% limite of results to plot
imatch = (dataplot.ninsert<=o.nmolplot);
%%%% plot concentration results 
for i = 1:size(dataplot.Caveragetab,1)
        x = dataplot.xplot(imatch);
        y = dataplot.Ctabout(imatch,i);
        y(y<o.ylim) = o.ylim;
        text(x,y,marker{i},'fontsize',o.markersize,'color',colmol(i,:),'VerticalAlignment','middle','HorizontalAlignment','center')
end
    
if o.plotfit          
    for i = 1:o.nplotfit     
        x = dataplot.xplot(imatch);
        y = dataplot.Ctabout(imatch,i);
        y(y<o.ylim) = o.ylim ;
        isegment = 1:length(y); isegment(y<o.ylim)=Inf;
        segment = grpseq(isegment);
        % fit C by segments
        if size(segment,1) > 1
            for ifit = 1:length(segment)-1;
                valid = ~isnan(y(segment{ifit}));
                if length(y(~isnan(y(segment{ifit})))) > 1
                    ylog = log(y(segment{ifit}(valid)));
                    sp = spaps(x(segment{ifit}(valid)),ylog,0.01*var(ylog)*length(y(segment{ifit}(valid))));
                    xp = linspace(min(x(segment{ifit}(valid))),max(x(segment{ifit}(valid))),1000);
                    plot(xp,exp(fnval(sp,xp)),'color',colmol(i,:))
                else
                    plot(x(segment{ifit}),y(segment{ifit}),'color',colmol(i,:))
                end                   
            end   
        else
            xp = linspace(min(x(~isnan(y))),max(x(~isnan(y))),100);
            sp = csaps(x,y,[],xp);
            plot(xp,sp,'color',colmol(i,:))
        end
    end
end
% format
xlabel('Number of inserted substances','fontsize',o.fontsize1)
ylabel(o.ylabel,'fontsize',o.fontsize1)
formatax(hs(1),'xlim',[0 o.nmolplot+1],'ylim',[o.ylim dbout.valuesmax*2],'fontsize',o.fontsize1,'yscale',o.yscale,'linewidth',1)

imolmax = max(max(dbout.indmol{o.nmolplot}));
%%%% get add data for likely conc
nscenario = zeros(imolmax,1);
likely = cell(imolmax,1);
likelyconc = zeros(imolmax,1);
likelyconcCI = zeros(imolmax,1);
for i = 1:imolmax
    y = dataplot.Ctabout(imatch,i);    
    y(y<o.ylim) = o.ylim ;
    % get number of scenarios
    nscenario(i) = length(dataplot.Ctabout(~isnan(y),i));
    % get likely value (mean value of longest segment)
    isegment = 1:length(y); isegment(y<o.ylim)=Inf;
    segment = grpseq(isegment);
    if size(segment,1) > 1
        for ifit = 1:length(segment)-1;
            valid = ~isnan(y(segment{ifit}));
            likely{i}.lenght(ifit) = length(y(segment{ifit}(valid)));
            likely{i}.value(ifit) = mean(y(segment{ifit}(valid)));
            likely{i}.CI(ifit) = 1.96*std(y(segment{ifit}(valid)))/sqrt(length(y(segment{ifit}(valid))));
        end 
    else
        valid = ~isnan(y(segment{1}));
        likely{i}.lenght = length(y(segment{1}(valid)));
        likely{i}.value = mean(y(segment{1}(valid)));
        likely{i}.CI = 1.96*std(y(segment{1}(valid)))/sqrt(length(y(segment{1}(valid))));
    end
    [~,imax] = max(likely{i}.lenght);
    likelyconc(i) = ceil(likely{i}.value(imax)); 
    likelyconcCI(i) = ceil(likely{i}.CI(imax));
end
    
%%%%% add legend
subplot(hs(2)), hsleg = subplots(1,[.5 .5],0,.025,'strict');
subplot(hsleg(1)), hs2 = subplots([.1 .4 .1 .15 .25],[.38 .62],0,0);
subplot(hs2(1,1)), text(.5,.5,'Class','fontsize',o.fontsize2,'VerticalAlignment','middle','HorizontalAlignment','center')
subplot(hs2(1,2)), text(.5,.5,sprintf('%s\n%s','Likely','substances'),'fontsize',o.fontsize2,'VerticalAlignment','middle','HorizontalAlignment','center')
subplot(hs2(1,3)), text(.5,.5,'\rho','fontsize',o.fontsize2,'VerticalAlignment','middle','HorizontalAlignment','center')
subplot(hs2(1,4)), text(.5,.5,sprintf('%s\n%s','Number of','scenarios'),'fontsize',o.fontsize2,'VerticalAlignment','middle','HorizontalAlignment','center','rotation',90)
subplot(hs2(1,5)), text(.5,.5,sprintf('%s\n%s','Likely conc.','(mg\cdotkg^{-1})'),'fontsize',o.fontsize2,'VerticalAlignment','middle','HorizontalAlignment','center')

formatax(hs2(1,:),'ylim',[0 1],'xlim',[0 1],'xtick',[],'xticklabel','','ytick',[],'yticklabel','','linewidth',.75)

subplot(hs2(2,1))
for i = 1:imolmax
    text(0.5,imolmax-(i-1),sprintf('C%0.2d',o.class(i)),'color',colmol(i,:),'fontsize',o.markersize,'VerticalAlignment','middle','HorizontalAlignment','center')
end

subplot(hs2(2,2))
for i = 1:imolmax
    text(0.02,imolmax-(i-1),sprintf('%s: %s',marker{i},dataplot.substancename{i}),...
        'color',colmol(i,:),'fontsize',o.markersize,'VerticalAlignment','middle','HorizontalAlignment','left')
end

subplot(hs2(2,3))
for i = 1:imolmax
    text(.5,imolmax-(i-1),sprintf('%0.2g',o.rhosort(i)),'color',colmol(i,:),'fontsize',o.markersize,'VerticalAlignment','middle','HorizontalAlignment','center')
end
subplot(hs2(2,4))
for i = 1:imolmax
    text(.5,imolmax-(i-1),sprintf('%0.2g',nscenario(i)),'color',colmol(i,:),'fontsize',o.markersize,'VerticalAlignment','middle','HorizontalAlignment','center')
end

subplot(hs2(2,5))
for i = 1:imolmax
    if isnan(likelyconc(i))
        text(.5,imolmax-(i-1),'<DL','color',colmol(i,:),'fontsize',o.markersize,'VerticalAlignment','middle','HorizontalAlignment','center')
    else
        text(.5,imolmax-(i-1),sprintf('\\fontsize{8}%0.2g \\fontsize{6}\\pm%0.2g',likelyconc(i),likelyconcCI(i)),...
            'color',colmol(i,:),'VerticalAlignment','middle','HorizontalAlignment','center')    
    end
end

for i = 1:size(hs2,2)
    subplot(hs2(2,i)), line([0 1],[1 1]*(imolmax-ifalsepos+.5),'color','k','linewidth',1,'linestyle','--')
end
  
formatax(hs2(2,:),'ylim',[0.5 imolmax+.5],'xlim',[0 1],'xtick',[],'xticklabel','','ytick',[],'yticklabel','','linewidth',0.75) 
%%%% add combinationt tree
subplot(hsleg(2)), hs3 = subplots(1,[.9 .1],0,.02,'strict');
subplot(hs3(1)), hs3a = subplots([margingraph 1-2*margingraph margingraph],[.05 .85 .1],0,0,'alive',5,'strict');
subplot(hs3a), scfd(dataplot.graph.figuredata,'noaxes','nolegend'), axis tight, set(hs3a,'visible','off')
% colorbar
subplot(hs3(2)), hs3b = subplots([o.margincolorbar 1-2*o.margincolorbar o.margincolorbar],1,0,0,'alive',2); % note that color @jet should match the one used in nmrlsqnonneg (default = jet(64))
subplot(hs3b)
poscolorbar = get(gca,'position'); poscolorbar = [poscolorbar(1) poscolorbar(2)+.01 poscolorbar(3) poscolorbar(4)/4.5];
colorlim = (ceil(linspace(dbout.valuesmin,dbout.valuesmax,5)));
colorbarcustom('colormap',dataplot.graph.colormap,'position',poscolorbar,'fontsize',o.fontsize2,'LineWidth',.65,...
           'ctick',colorlim,'xlim',[colorlim(1) colorlim(end)],'xtick',colorlim,'xticklabel',colorlim)
xlabel('conc. (mg\cdotkg^{-1})','fontsize',o.fontsize1)       

set(hs3b,'visible','off')  