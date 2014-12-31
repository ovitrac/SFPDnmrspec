function hfigout = plotdeformulateNMR(R,varargin)
%PLOTDEFORMULATENMR plot results of deformulation by NMR technique using output of deformulateNMR or nmrpairidentification, nmrdeconvolution
% To plot concentration (plotnmrlsqnonneg): use output of deformulateNMR or nmrdeconvolution
% To plot pairwise correlation coefficient: use output of deformulateNMR or nmrpairidentification
%
% SYNTAX
% plotdeformulateNMR(varargin) 
% VARARGIN
%          'figname': name of figure for save (default='mixturename_NMR_deconvolution')
%    'paperposition': paper position for figure (default=[3.0301    1.2380   23.6173   18.5080])
% 'paperorientation': paperorientation (default='landscape')
%     'marginlayout': margin to construct layout of report (default=[.2 .1 .2]
%     'spacerheater': spacer between areas to plot of heater layout (default=0.05);
%     'spacerfooter': spacer between areas to plot of footer layout (default=0.08);
%   'margincolorbar': numeric (<1) for set axis subplot to plot colorbar (default = 0.1)
%      'margingraph': numeric (<1) for set axis subplot to plot tree graph of scenarios (default = 0.1)
%       'headertext': additional text in header area of layout (default=empty)
%       'footertext': additional text in footer area of layout (default=empty)
%       'resolution': resolution for saving figure (default=300)
%       'outputpath': string, directory of output for saving figures
% See also: nmrpairidentification, nmrdeconvolution, deformulateNMR, plotnmrlsqnonneg

% RMNSPEC v 0.5 - 30/09/2014 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 28/11/2014
% History
% 15/10/14: use localname only if outputpathdefault is empty
% 25/11/14: use uistack to move figure to the top 
%             remove 'Warning: Unable to interpret TeX string' for footer text and header text
% 28/11/14: add output for figure handle
% default
outputpathdefault = '';
default = struct('figname','',...
                 'paperposition',[3.0301    1.2380   23.6173   18.5080],...
                 'paperorientation','landscape',...
                 'marginlayout',[.2 .1 .2],...
                 'spacerheater',0.05,...
                 'spacerfooter',0.08,...
                 'margincolorbar',0.1,...
                 'margingraph',[],...
                 'headertext','',...
                 'footertext','',...
                 'fontsizeheadertext',12,...
                 'fontsizefootertext',14,...
                 'resolution',300,...
                 'outputpath',outputpathdefault);

default.help.figname = 'name of figure for save (default=mixturename_NMR_deconvolution)';
default.help.paperposition = 'paper position for figure (default=[3.0301    1.2380   23.6173   18.5080])';
default.help.paperorientation = 'paperorientation (default=landscape)';
default.help.marginlayout = 'margin to construct layout of report (default=[.2 .1 .2]';
default.help.spacer = 'spacer between areas to plot of layout (default=0.08)';
default.help.margincolorbar = 'numeric (<1) for set axis subplot to plot colorbar (default = 0.1)';
default.help.margingraph = 'numeric (<1) for set axis subplot to plot tree graph of scenarios (default = 0.1)';
default.help.resolution = 'resolution for saving figure (default=300)';
default.help.headertext = 'additional text in header area of layout (default=empty)';
default.help.footertext = 'additional text in footer area of layout (default=empty)';
default.help.outputpath = 'string, directory of output for saving figures';

keyword = {'plotconc', 'plotrho', 'print'};

% argcheck
o = argcheck(varargin,default,keyword,'case');
if nargin < 1, error('1 argument is required. Please use functions NMRPAIRIDENTIFICATION, NMRDECONVOLUTION or DEFORMULATENMR to get result'), end
if ischar(R)
    if exist(R,'file')
        dispf('\tload data to plot'), fileinfo(R), load(R)
    else
        error('the data to plot ''%s'' does not exist in ''%s''',lastdir(R),rootdir(R))
    end
end
if ~isempty(o.figname), figname = o.figname; 
elseif o.plotconc
    figname = sprintf('%s_NMR_deconvolution',R.mixturename); 
elseif o.plotrho
    figname = sprintf('%s_NMR_rho_ranking',R.mixturename);
end

if ~isempty(o.headertext), headertext = regexprep(o.headertext,{'\' '_'},{'\\\' '\\_'});
elseif o.plotconc
    headertext = sprintf('\\fontsize{14}Sample name: %s\n\\fontsize{11}',regexprep(R.mixturename,{'\' '_'},{'\\\' '\\_'}));
elseif o.plotrho
    headertext = sprintf('\\fontsize{14}Sample name: %s\n\\fontsize{11}',regexprep(R.mixturename,{'\' '_'},{'\\\' '\\_'}));
end
if ~isempty(o.footertext), footertext = regexprep(o.footertext,{'\' '_'},{'\\\' '\\_'});
else footertext = sprintf('Database: %s \nDate: %s',regexprep(fullfile(R.dbinfo.path,R.dbinfo.filename),{'\' '_'},{'\\\' '\\_'}),regexprep(R.dbinfo.date,{'\' '_'},{'\\\' '\\_'}));
end
if isempty(o.outputpath)
    switch localname
        case 'LP-MOL5'
            o.outputpath = 'C:\Data\Olivier\INRA\Projects\SafeFoodPack_design\livrables\D1_deformulation_rule';
        case 'LP_MOL2'
            o.outputpath = 'C:\Data\Olivier\INRA\Projects\SafeFoodPack_design\livrables\D1_deformulation_rule';          
    end
end
outputpath = fullfile(o.outputpath,R.mixturename);
if ~exist(outputpath,'dir'), mkdir(outputpath); end
% close all
% plot
hout = SFPDlayout('figname',figname,'paperposition',o.paperposition,'paperorientation',o.paperorientation,'margin',o.marginlayout,'spacerheater',o.spacerheater ,'spacerfooter',o.spacerfooter);
if o.plotconc % plot concentration
    subplot(hout.center), plotnmrlsqnonneg(R.dataplot,R.dbout,'parent',hout,'class',R.rang,'rhosort',R.rhotestsort,'plotfit','nplotfit',R.nplotfit,'margincolorbar',o.margincolorbar,'margingraph',o.margingraph)
    subplot(hout.header), text(.5,.5,headertext,'fontsize',o.fontsizeheadertext,'verticalalignment','middle','horizontalalignment','center')
    subplot(hout.footer), text(.01,.5,footertext,'fontsize',o.fontsizefootertext,'verticalalignment','middle','horizontalalignment','left')
elseif o.plotrho % plot ranking of pairwise correlation coefficient
    nmol = length(R.corrmeanweightsort);
    subplot(hout.center)
    hstmp = subplots([.2 .8],1,0,0,'alive',2);
    subplot(hstmp(1)), plot(R.corrmeanweightsort,1:nmol,'-o','linewidth',1)
    xlabel('\rho^2','fontsize',14)
    formatax(gca,'ylim',[0 nmol+1],'ytick',0:nmol+1,'xlim',[0 1],'xtick',0:0.1:1,'yticklabel','','fontsize',14)
    for i = 1:nmol
        text(-.01,i,R.molsortbyrho{i},'fontsize',8,'VerticalAlignment','middle','HorizontalAlignment','right')
    end
    subplot(hout.header), text(.5,.5,headertext,'fontsize',o.fontsizeheadertext,'verticalalignment','middle','horizontalalignment','center')
    subplot(hout.footer), text(.01,.5,footertext,'fontsize',o.fontsizefootertext,'verticalalignment','middle','horizontalalignment','left')
end
% save
% saveas(gcf,fullfile(outputpath,[figname '.fig']))
if o.print
    uistack(hout.hfig,'top')
    print_pdf(o.resolution,figname,outputpath,'nocheck','PaperPositionMode','auto')
    print_png(o.resolution,figname,outputpath,'',0,0,0)
end

if nargout, hfigout = hout.hfig; end