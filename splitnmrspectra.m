function [x,y] = splitnmrspectra(mollist)
%SPLITNMRSPECTRA splits NMR spectrum into several sub-spectra 
% Syntax: y = splitnmrspectra(mollits)
%
% Example
% mol = {'Erucamide','Tinuvin326'};
% [x,y] =  splitnmrspectra(mol);
% for i = 1:length(mol)
%     figure, plot(x,y(:,:,i)), title(mol(i))
% end

% RMNSPEC v 0.1 - 18/03/13 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 
% history

% argcheck
if ~iscell(mollist), mollist = {mollist}; end
if ~iscellstr(mollist), error('mollist must be a string'), end
% main
% load nmr databases
[dbpur,dbxpur] = nmrloadbase;

subdb = nmrsubdb(dbxpur,'commonname',mollist);  % extract sub-base of spectra 
nmol = length(mollist);
x = dbxpur.ppm;
for imol = 1:nmol
    ROI = dbpur.(mollist{imol}).gates; % extract RoI (nx4 array: ppmmin, ppmmax, buffer and weight)
    nROI = length(ROI);
%     y = cell(dbxpur.m,nROI,nmol);
    for i = 1:nROI
        novalid = ((x<=ROI(i,1))|(x>=ROI(i,2))); % other areas than considered ROI
        tmp = subdb.I0(:,imol);
        tmp(novalid) = 0; % 0 for other areas than considered ROI
        y(:,i,imol) = tmp*ROI(i,4); % y(:,i,imol) = {tmp(:,1)*ROI(i,4)};  multiplying of weight to remove excluded areas (in mask database, exluded areas (standards, parasite peaks...) are weighted by 0)
    end
end

% default
% default = struct(...
%                 'dbxpur',[],...
%                 'dbmask',[]...
%                 );
% argcheck
% o = argcheck(varargin,default);
% if ~iscell(mollist), mollist = {mollist}; end
% if ~iscellstr(mollist), error('mollist must be a string'), end
% if isempty(o.dbxpur) && evalin('base','exist(''dbxpur'',''var'')')
%    o.dbxpur = evalin('base','dbxpur');
% end
% if isempty(o.dbmask)
%     if evalin('base','exist(''dbmask'',''var'')')
%         o.dbmask = evalin('base','dbmask');
%     end
% end