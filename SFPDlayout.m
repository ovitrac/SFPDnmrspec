function h = SFPDlayout(varargin)
%SFPDlayout returns handles for layout of report for SFPD project with logo and others information
% SYNTAX
% h = SFPDlayout()
% [h,hs] = SFPDlayout()
% VARARGIN
%              'figname': figure's name (default = 'report')
%        'paperposition': position of figure, (default = [])
%            'papertype': type of paper (A3,A4...), (default = 'A4')
%     'paperorientation': paper orientation (portrait or landscape) (default='portrait')
%               'margin': 1x3 vector [header footer leftcorner] for (default = [.12 .06 .2])
%         'spacerheater': gap/space between heater and main area to plot
%         'spacerfooter': gap/space between footer and main area to plot
% OUTPUT : structure h with fields
%                   hfig: figure handle;
%                 center: handle of principal axe subplot
%                 header: handle of axe for header
%                 footer: handle of axe for footer

% rmnspec v
% 15/10/2014 - INRA\Olivier Vitrac, Mai Nguyen - rev. 20/11/14
% history
% 20/11/14: change output s structure

% default
default = struct('figname','report','paperposition',[],'paperorientation','portrait','papertype','A4',...
                  'margin',[.12 .06 .2],'spacerheater',0.05,'spacerfooter',0.05);

% argcheck
o = argcheck(varargin,default); 

% general axes for plot
figname = o.figname;
hfig = figure; formatfig(hfig,'figname',figname,'paperposition',o.paperposition,'paperorientation',o.paperorientation,'papertype',o.papertype)
hstmp = [subplots([o.margin(3) 1-o.margin(3)],[o.margin(1) 1-(o.margin(1)+o.margin(2)) o.margin(2)],0,[o.spacerheater o.spacerfooter],'alive',[1 3 4 6]);
         subplots(1,[o.margin(1) 1-(o.margin(1)+o.margin(2)) o.margin(2)],0,[o.spacerheater o.spacerfooter],'alive',2,'position',gcf)];
set(hstmp,'xtick',[],'ytick',[])

% add logo and info for layout
currentfile = mfilename('fullpath');
currentpath = rootdir(currentfile);
subplot(hstmp(1)), hstmp1 = subplots([.02 .73 .25],1,0,0,'position',hstmp(1),'alive',2:3);
subplot(hstmp1(1)), imshow(imread(fullfile(currentpath,'private','logo_layout.png')))
subplot(hstmp1(2)), text(.5,.05,{'ANR-10-ALIA-009'},'verticalalignment','middle','horizontalalignment','left','fontsize',10,'rotation',90)
set(hstmp1,'xtick',[],'ytick',[],'visible','off')


info = sprintf('\\bfDate: \\rm%s\n\\bfAuthor: \\rm%s\n\\bfEngine: \\rm%s',date,username,regexprep(localname,'_','\\_'));
subplot(hstmp(2)), text(.01,0.5,info,'verticalalignment','middle','horizontalalignment','left','fontsize',10)    
    
% output
if nargout
    h.hfig = hfig;
    h.center = hstmp(5);
    h.header = hstmp(3);
    h.footer = hstmp(4);
end
