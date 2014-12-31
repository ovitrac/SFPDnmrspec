function subdbx = nmrsubdb(dbx,field,keys)
%NMRSUBDB extracts a subdatabase from a numerical database (dbx) loaded with NMRLOADASCII

% syntax: out = nmrsubdb(dbx,field,keys)
% dbx: database of NMR spectra (for pur subtances or mixtures) for mathematical treatment with fields (created by NMRLOADASCCI)
% output : subdbx contains same fields as dbx (input) created by NMRLOADASCII

% See also: nrlmoadascii

% example : use subdbxpur and subdbxmix to decompose a reference mixture (SFPDPP3)
% ----------------------------------------------------------
%   [dbpur,dbxpur] = nmrloadascii('odsfile','substance.ods','sheetname','substance','path','data_pur','ppmmin',-0.5,'ppmmax',12);
%   [dbmix,dbxmix] = nmrloadascii('odsfile','extract.ods','sheetname','extract','primarykey','reference','path','data_mixture','ppmmin',-0.5,'ppmmax',12);
%   % create a subdbxpur
%   substoextract = {'MBOCA156',...
%                  'MBOCA313',...
%                  'MBOCA625',...
%                  'MBOCA1250',...
%                  'MBOCA3125'...
%                  };
%   subdbxMBOCA = nmrsubdb(dbxpur,'commonname',substoextract);  
% 
%   %  nmrcorr for pair of substances
%   corr = nmrcorr(subdbxMBOCA,subdbxMBOCA,'all');
%    
%   % plot 
%   for n = 1:length(corr.rho)  
%       [~,imax] = max(corr.Cij(:,n));
%       figure
%       hs = subplots([0.65 0.35],1,0.05,0);     
%       subplot(hs(1)), plot(subdbxMBOCA.ppm+corr.lag(imax),subdbxMBOCA.I(:,corr.i(n)),subdbxMBOCA.ppm,subdbxMBOCA.I(:,corr.j(n))), xlim([-0.5 12])
%       ylabel('Intensity','fontsize',14)
%       set(gca,'XDir','reverse')
%       xlabel('chemical shift (ppm)','fontsize',14)    
%       legend(subdbxMBOCA.commonname([corr.i(n) corr.j(n)]),'location','NorthWest')
%       title(['correlation coefficient: ' num2str(corr.rho(n))])
%       plotdata = gcfd;
%       % add zoom
%       hzoom = [axes('position',[0.3514    0.6597    0.1712    0.2433])
%                axes('position',[0.1430    0.4317    0.1791    0.2510])];
%       subplot(hzoom(1)), scfd(plotdata,'noaxes','nolegend'); axis([2  4.2   -0.5e-8   5e-8]), set(gca,'XDir','reverse')
%       subplot(hzoom(2)), scfd(plotdata,'noaxes','nolegend'); axis([6.5  8.1   -0.5e-8   5e-8]), set(gca,'XDir','reverse')
%       set(hzoom,'box','on','fontsize',6)
%     
%       subplot(hs(2)), plot(corr.lag(:,1),corr.Cij(:,n))
%       figname = [subdbxMBOCA.commonname{corr.i(n)} '_' subdbxMBOCA.commonname{corr.j(n)}];
%       formatfig(gcf,'paperposition',[0.6694    1.7960   28.3387   17.3920],'PaperOrientation','landscape')
%       print_pdf(300,figname,figurefolder,'nocheck','PaperPositionMode','auto')
%       print_png(300,figname,figurefolder)
%   end
% ---------------------------------------------------------------------------------------------------------------------------


% RMNSPEC v 0.1 - 10/09/12 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 11/09/12

% History
% 11/09/12 : add help and example

% arg check
if nargin<3, error('3 arguments are required'), end
if ~isstruct(dbx) || ~isfield(dbx,'I') || ~isfield(dbx,'ppm') || ~isfield(dbx,'m') || ~isfield(dbx,'n') || ~isfield(dbx,'step')
    error('dbx must be created by NMRLOADASCII')
end
if ~isfield(dbx,field), error('the field ''%s'' does not exist in the supplied database',field), end
if ~ischar(dbx.(field)) && ~iscellstr(dbx.(field)), error('the field ''%s'' in the supplied database must contains only strings ',field), end
if ~iscell(keys), keys = {keys}; end
if ~iscellstr(keys), error('the keys must be given in a cell array of strings (numerical values are excluded)'), end
nkeys = length(keys);
if nkeys>length(unique(keys)), error('keys must be unique, please check'), end
[~,ifield] = intersect(dbx.(field),keys);
if length(ifield)<nkeys
    dispf('ERROR\t%d values given as keys cannot be found in ''%s''',nkeys-length(ifield),field)
    cellfun(@(m) dispf('\t''%s'' is missing',m),setdiff(keys,dbx.(field)))
    error('%d missing values (see above)',nkeys-length(ifield))
end

% extract database
subdbx = dbx;
subdbx.n = nkeys;
for f = fieldnames(dbx)'
    if isnumeric(dbx.(f{1})) && (size(dbx.(f{1}),2)==dbx.n)
        subdbx.(f{1}) = dbx.(f{1})(:,ifield);
    elseif iscell(dbx.(f{1})) && (numel(dbx.(f{1}))==dbx.n)
        subdbx.(f{1}) = dbx.(f{1})(ifield);
    end
end
