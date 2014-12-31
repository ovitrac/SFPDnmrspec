function W = multipb(x,gates)
% MULTIPB builds n gates using pb: W = sum(i=1...N)(weight(i)*pb(x,[xstart xfinal],s))
% SYNTAX 
%   p = multipb(x,gates)
% INPUTS
%      x: mx1 array (scale) or dbx spectral databse created with NMRLOADASCII
%  gates: nx4, nx3 or nx2 array with
%         column 1: starting x value
%         column 2: final x value
%         column 3: buffer (same units as x)
%         column 4: corresponding weight 
%
% OUTPUT
%      W: mx1 array, gates weighted
%
% See also: pb
%
% Example
% [dbpur,dbxpur] = nmrloadascii('path', fullfile(find_path_toolbox('rmnspec'),'data_pur'),'ppmmax',13);
% dbmask = nmrloadmask('path', fullfile(find_path_toolbox('rmnspec'),'data_pur'));
% gates = dbmask.Tinuvin326.gates;
% W = multipb(dbxpur.ppm,gates);
% figure, hs = subplots(1,[1 1 1],0,0);
% subplot(hs(1)), plot(dbxpur.ppm,dbxpur.I(:,27)), xlim([-1 12]), ylabel('Intensity','fontsize',14)
% title('Tinuvin 326')
% subplot(hs(2)), plot(dbxpur.ppm,W), xlim([-1 12]), ylabel('Weight','fontsize',14)
% subplot(hs(3)), plot(dbxpur.ppm,dbxpur.I(:,27).*W), xlim([-1 12]), ylabel('Intensity','fontsize',14)
% set(hs,'XDir','reverse')   
% xlabel('chemical shift (ppm)','fontsize',14)
% 
% figure, hs = subplots(1,[1 1],0,0);
% subplot(hs(1)), plot(dbxpur.ppm,dbxpur.I(:,27)); 
% ax1 = gca;
% formatax(ax1,'box','off','ylim',[-.05e-7 max(dbxpur.I(:,27))],'xlim',[-1 12],'fontsize',10,'Xcolor','b','Ycolor','b','LineWidth',1)
% ylabel(ax1,'Intensity','fontsize',12)    
% ax2 = axes('position',get(ax1,'position'),'xAxisLocation','top','yAxisLocation','right','color','none','XColor','k','YColor','k','XDir','reverse'); %#ok<LAXES>
% hold on
% plot(dbxpur.ppm,W,'k','parent',ax2)
% formatax(ax2,'box','off','xlim',[-1 12],'ylim',[-.05 ceil(max(W))],'fontsize',10,'LineWidth',1), ylabel(ax2,'Weight','fontsize',12)    
% title('Tinuvin 326','fontsize',14)
% subplot(hs(2)), plot(dbxpur.ppm,dbxpur.I(:,27).*W), 
% formatax(hs(2),'xlim',[-1 12],'fontsize',10,'LineWidth',1), ylabel('Intensity','fontsize',12)
% xlabel('chemical shift (ppm)','fontsize',12)
% set(hs,'XDir','reverse')  
% datacopy = gcfd;
% 
% hzoom = [axes('position',[0.4199    0.6828    0.3074    0.2046],'color','none')
%          axes('position',[0.4199    0.6828    0.3074    0.2046],'color','none','yAxisLocation','right')
%          axes('position',[0.4189    0.2678    0.3074    0.2046])];
% hold on
% subplot(hzoom(1)), scfd(datacopy(1,1),'noaxes','nolegend'); axis([7.1  8.2   0   1e-7])
% subplot(hzoom(2)), scfd(datacopy(1,3),'noaxes','nolegend'); axis([7.1  8.2   0   4])
% subplot(hzoom(3)), scfd(datacopy(1,2),'noaxes','nolegend'); axis([7.82  8.12   -.05e-7   1.5e-7])
% set(hzoom,'box','on','XDir','reverse')
% 
% formatfig(gcf,'paperposition',[0.7530    0.7975   28.1714   19.3890],'PaperOrientation','landscape')
% figname = 'tinuvin326';
% print_pdf(300,figname,fullfile(find_path_toolbox('rmnspec')),'nocheck','PaperPositionMode','auto')
% print_png(300,figname,fullfile(find_path_toolbox('rmnspec')))
%
% Advanced example
% x = 0:100; 
% gates = [10 12; 25 29; 50 55];
% W = multipb(x,gates);
% figure, plot(x,W)
%
% RMNSPEC v 0.1 - 23/01/13 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev.
% 
% History 
%
% argcheck
if isstruct (x) && isfield(x,'ppm') && isfield(x,'I') && isfield(x,'m') && isfield(x,'n')
    if nargin<2, gates = x.gates; end
    x = x.ppm;
elseif nargin<2, error('2 arguments are required')
end
x = x(:);
m = length(x);
n = size(gates,1);
% main
if size(gates,2)<2, error('gates should be at least a nx2 array')
elseif size(gates,2)>4, error('gates should be at maximum a nx4 array'), end
if size(gates,2)<3, gates = [gates ones(n,1)/100]; end
if size(gates,2)<4, gates =[gates ones(n,1)]; end
P = zeros(m,n);
W = zeros(m,1);
for i = 1:n
    P(:,i) = pb(x,[gates(i,1) gates(i,2)],[1 -1]*gates(i,3));
end
W = P*gates(:,4);

