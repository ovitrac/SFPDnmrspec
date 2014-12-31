function [spectra,ysubs] = nmrbuildtheospec(listsubstances,dbfit,dbxpur,weight,varargin)
%NMRBUILDTHEOSPEC builds/ generates theoretical spectrum of mixture of fitted/theoretical substances
% syntax: spectra = nmrbuildtheospec(listsubstances,dbfit,[],weight)
% INPUTS            
%   listsubstances: 1 x n cell {'mol1' 'mol2' 'mol3'}, n substances to build mixture (theoretical, fitted substances)
%            dbfit: theoretical , fitted NMR database containing substances
%           dbxpur: nmr database containing real substance spectra data
%           weight: m x 1 array, [w1;w2;...;wm] weight of m bands (motifs) of all subtances in mixture
%                   n x 1 cell {[w11;w21;w31];[w12;w22;w32];...} weight of all bands (motifs) of n subtances in mixture stored by susbtance
%                   n x 1 array, [w1;w2;...wn] weight of n substances (all bands in the same substance have same weight)
% VARARGIN
%        dbfitpath: location of fitting NMR spectra database (default = fullfile(find_path_toolbox('rmnspec'),'data_pur'))
%        dbfitname: fitting NMR spectra database (default = dbfit.mat, 
%          ppmarea: [ppm1 ppm2 k] array, gate to extract band of spectra
%                   ppm1: inferior borne (default=-1)
%                   ppm2: superior borne (default=13)
%                      k: number of points in gate (default=2^15)
% KEYWORD
%           'Ifit': theoretical signals to be use
%           'Iraw': real signals to be used
%     'normalized': normalization of signals by number of H and pente(calibration data)

% OUTPUTS
% spectra: resolution x 2 array, spectrum of mixture built
%          [:,1]: xscale, ppm
%          [:,2]: yscale, intensity
%
% EXAMPLE 1
% dbfit = nmrloaddbfit('path',fullfile(find_path_toolbox('rmnspec'),'data_pur'),'dbname','dbfit.mat');
% [dbpur,dbxpur] = nmrloadbase;
% mol = {'Erucamide' 'Tinuvin326' 'Irganox1076'}; 
% narea = 0; for i=1:length(mol), narea = narea + size(dbfit.(mol{i}).all,1); end
% weight = zeros(narea,1); k = 0; 
% for i=1:length(mol), 
%     ROI = dbpur.(mol{i}).gates; ROI = ROI(ROI(:,4)>0,:);
%     for j = 1:length(ROI), k=k+1; weight(k,1) = ROI(j,4); end
% end
% mixsp = nmrbuildtheospec(mol,dbfit,dbxpur,weight); 
%
% EXAMPLE 2
% dbfit = nmrloaddbfit('path',fullfile(find_path_toolbox('rmnspec'),'data_pur'),'dbname','dbfit.mat');
% [dbpur,dbxpur] = nmrloadbase;
% mol = {'Erucamide' 'Tinuvin326' 'Irganox1076'}; 
% weight = cell(length(mol),1);
% for i=1:length(mol) 
%     ROI = dbpur.(mol{i}).gates; ROI = ROI(ROI(:,4)>0,:);
%     weight{i,1} = ROI(:,4);
% end
% mixsp = nmrbuildtheospec(mol,dbfit,dbxpur,weight,'ppmarea',[4 9 5000]); 
%
% EXAMPLE 3
% dbfit = nmrloaddbfit('path',fullfile(find_path_toolbox('rmnspec'),'data_pur'),'dbname','dbfit.mat');
% [dbpur,dbxpur,dbmask,mask,dbmultiplet,dbmix,dbxmix] = nmrloadbase;
% mol = {'Oleamide' 'Erucamide' 'Chimassorb944' 'Tinuvin622'};
% weight = [1 2 4 6];
% [mixsp,ysubstance] = nmrbuildtheospec(mol,dbfit,dbxpur,weight); % theoretical mixture
% dbcorrsub = nmrloadmatrix(mol,[],dbxpur,dbfit); % matrix of correlation between substances
% dbcorrmix = nmrloadmatrix(mol,mixsp(:,2),dbxpur,dbfit); % matrix of correlation between substances and theoretical mixture
% W = lsqlin(dbcorrsub.corrzerolag,dbcorrmix.corrzerolag,[],[],[],[],0);
% [~,indmol] = intersect(mol,dbcorrsub.subtestlist);
% Itoplot = cellfun(@(m) dbpur.(m).In,dbcorrsub.subtestlist,'UniformOutput',false);
% Itoplot = cat(2,Itoplot{:});
% figure, plot(mixsp(:,1),[mixsp(:,2) ysubstance(:,indmol)*W]), legend({'reference mixture' 'deconvolution solutions'}) % Itoplot*W
% for i = 1:length(mol)
%     figure, plot(mixsp(:,1),[ysubstance(:,indmol(i))*weight(indmol(i)) ysubstance(:,indmol(i))*W(i)]) % Itoplot(:,i)*W(i)
%     legend({'reference' 'deconvolution solutions'}) 
%     title(dbcorrsub.subtestlist{i})
% end
% See also: nmrloadmatrix
%
% RMNSPEC v 0.1 - 16/07/2013 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev.22/07/13
%
% History
% 22/07/13 add example
% 19/10/13 add keyword Ifit' or Iraw
% default
default = struct('dbfitpath',fullfile(find_path_toolbox('rmnspec'),'data_pur'),...
                 'dbfitname','dbfit.mat',...
                 'ppmarea',[-1 13 2^15],...
                 'dbcalibration',[]); 
keyword = {'Ifit','Iraw','normalized'}; 

% arg check
o = argcheck(varargin,default,keyword,'nostructexpand');
if nargin<1, error('1 argument is required'), end
if nargin<2, dbfit = []; end
if nargin<3, dbxpur = []; end
if nargin<4, weight = []; end
if isempty(dbfit), dbfit = nmrloaddbfit('path',o.dbfitpath,'dbname',o.dbfitname); end
if isempty(dbxpur), [~,dbxpur] = nmrloadbase; end % load dbxpur NMR spectra database
if o.normalized && isempty(o.dbcalibration)
    o.dbcalibration = nmrloaddbspec('dbcalibration');
    nfo = fileinfo(fullfile(find_path_toolbox('rmnspec'),'nmrdatabase.mat'));
    fprintf('DBCALIBRATION is not defined, use dbcalibration in following base:\n%s\ndate and bytes:%s%s',strcat(nfo.path,nfo.filename),nfo.date,nfo.bytes)
end
if iscell(weight), weight = cell2mat(weight(:)); end % convert weight in array
if isvector(weight), weight = weight(:); end
if size(weight,1) < length(listsubstances), error('Check size of weight! Weight must be defined in array or cell'), end
narea = 0; for i=1:length(listsubstances), narea = narea + size(dbfit.(listsubstances{i}).all,1); end
if  size(weight,1) > length(listsubstances) &&  size(weight,1) < narea, error('Check size of weight! Weight must be defined in array or cell'), end
ppm = dbxpur.ppm; m = dbxpur.m; step = dbxpur.step;

commonname = fieldnames(rmfield(dbfit,'help'));
if o.Ifit, signaltype = 'Ifit'; end
if o.Iraw, signaltype = 'Iraw'; end

% check list of substances
if ~iscell(listsubstances), listsubstances = {listsubstances}; end
if ~iscellstr(listsubstances), error('the substance list must be given in a cell array of strings (numerical values are excluded)'), end 
nlistsubstances = length(listsubstances); 
if nlistsubstances>length(unique(listsubstances)), error('substance names must be unique, please check'), end
[subcommonname,~] = intersect(commonname,listsubstances); % name of substances in alphabetic order =/ than order in dbfit
if length(subcommonname)<nlistsubstances
    dispf('ERROR\t%d values given as substance names cannot be found',nlistsubstances-length(subcommonname))
    cellfun(@(m) dispf('\t''%s'' is missing',m),setdiff(listsubstances,commonname))
    error('%d missing values (see above)',nlistsubstances-length(subcommonname))
end

% main
% -> extract substance spectra for building mixture spectrum
% -> use original listsubstances to keep original order
I = zeros(m,nlistsubstances); % 
k=0;
for i = 1:nlistsubstances
    nROI = size(dbfit.(listsubstances{i}).all,1);
    ywhole = zeros(m,nROI);   
    for iROI = 1:nROI
        k=k+1;
        if o.normalized, y = (dbfit.(listsubstances{i}).all(iROI,1).(signaltype))/(o.dbcalibration.dbcalibration.(listsubstances{i}).p(1)*sum(o.dbcalibration.dbcalibration.(listsubstances{i}).H));
        else y = (dbfit.(listsubstances{i}).all(iROI,1).(signaltype));
        end
        x = dbfit.(listsubstances{i}).all(iROI,1).ppm;
        firstsegment = zeros(round((min(x)-min(ppm))/step),1);
        lastsegment = zeros(round((max(ppm)-max(x))/step),1);
        if ~isempty(weight)   % if weight is defined            
            if numel(weight)==nlistsubstances % weight is defined in array for each substance
                ywhole(:,iROI) = [firstsegment;y;lastsegment]*weight(i,1);
            else                              % weight is defined in array for each band
                ywhole(:,iROI) = [firstsegment;y;lastsegment]*weight(k,1);
            end
        else                   % no weight -> weight = 1;
            ywhole(:,iROI) = [firstsegment;y;lastsegment]; 
        end
    end
    I(:,i) = sum(ywhole,2);
end

% Output
spectra = zeros(o.ppmarea(3),2);
tmpsp = sum(I,2);
spectra(:,1) = linspace(o.ppmarea(1),o.ppmarea(2),o.ppmarea(3))';
spectra(:,2) = interp1(ppm,tmpsp,spectra(:,1));
if nargout >1, ysubs = I; end

%% plot 
if ~nargout
    yplot = zeros(m,nlistsubstances);
    for i= 1:nlistsubstances
        yplot(:,i) = I(:,i)+i*max(I(:))/10;
    end
    figure, plot(ppm,yplot)
    hold on, plot(spectra(:,1),spectra(:,2),'k')
    xlabel('chemical shift (ppm)','fontsize',12), ylabel('Intensity','fontsize',12)
    formatax(gca,'xlim',[o.ppmarea(1) o.ppmarea(2)]) 
    legend([listsubstances 'mixture spectrum'],'fontsize',10)
end