function [snr,mfiltout] = signaltonoise(X,mfilt)
%SIGNALTONOISE calculates the signal-to-noise ratio by comparing filtered (with filtzero) with non-filtered data
%  Syntax: snr = signaltonoise(X [,mfilt])
%          [snr,mfiltout] = signaltonoise(...)
%  INPUTS
%       X: mxn signal array
%   mfilt: 1xmf array of mfilt values to be used with filtzero (default = unique(max(2,round(logspace(-6,-2),7)*m))) );
%  OUTPUTS
%     snr: mfxn
%mfiltout: used mfilt values
%
% EXAMPLE:
%   if isempty(find_path_toolbox('rmnspec')), error('install first the toolbox rmnspec'), end
%   [dbpur,dbxpur] = nmrloadascii('path', fullfile(find_path_toolbox('rmnspec'),'data_pur'),'noprefetch');
%   X = cellfun(@(m) dbpur.(m).I,fieldnames(dbpur),'UniformOutput',false); X = cat(2,X{:});
%   mfilt = [2 5 10 15 20 25 30];
%   snr = signaltonoise(X,mfilt); [~,rank] = sort(snr(1,:),'descend');
%   set(figure,'DefaultAxesColorOrder',cbrewer('div','RdYlGn',size(snr,2))), plot(mfilt,snr(:,rank))
%   xlabel('m_{filt}','fontsize',12), ylabel('SNR (db)','fontsize',12)
%
% EXAMPLE 2 (based on  the example above)
%   mol = {'Stearicacid' 'Stearicacid400' 'Stearicacid100' 'MBOCA1250'};
%   X = cellfun(@(m) dbpur.(m).I,mol,'UniformOutput',false); X = cat(2,X{:});
%   snr = signaltonoise(X,mfilt); [~,rank] = sort(snr(1,:),'descend');
%   set(figure,'DefaultAxesColorOrder',cbrewer('div','RdYlGn',size(snr,2))), plot(mfilt,snr(:,rank))
%   xlabel('m_{filt}','fontsize',12), ylabel('SNR (db)','fontsize',12)
%
%   See for additional information on SNR: http://en.wikipedia.org/wiki/Signal-to-noise_ratio

% RMNSPEC v 0.1 - 06/10/12 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 07/10/12
% 
% History
% 07/10/12 add examples

% arg check
if nargin<1, error('one argument is required'); end
if nargin<2, mfilt = []; end
[m,n] = size(X);
if isempty(mfilt), mfilt = unique(max(2,round(logspace(-6,-2,7)*m))); end
mf = numel(mfilt);

% main
snr = zeros(mf,n);
screen = '';
for i=1:mf
    screen = dispb(screen,'SNR calculating %0.3g%%...',100*i/mf);
    Xf = filtzero(X,mfilt(i));
    snr(i,:) = 10*log10( var(Xf)./var(X-Xf) ); %10*log10( signal / noise );
end

snr(snr<0)=NaN;
if nargout>1, mfiltout = mfilt; end