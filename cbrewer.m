function [colormap]=cbrewer(ctype, cname, ncol, interp_method)
%
% CBREWER - This function produces a colorbrewer table (rgb data) for a 
% given type, name and number of colors of the colorbrewer tables. 
% For more information on 'colorbrewer', please visit
% http://colorbrewer2.org/
% 
% The tables were generated from an MS-Excel file provided on the website
% http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer_updates.html
%
% INPUT:
%   - ctype: type of color table 'seq' (sequential), 'div' (diverging), 'qual' (qualitative)
%   - cname: name of colortable. It changes depending on ctype.
%   - ncol:  number of color in the table. It changes according to ctype and
%            cname
%   - interp_method: if the table need to be interpolated, what method
%                    should be used for interp1? Default='cubic'.
% 
% A note on the number of colors: Based on the original data, there is
% only a certain number of colors available for each type and name of
% colortable. When 'ncol' is larger then the maximum number of colors
% originally given, an interpolation routine is called (interp1) to produce 
% the "extended" colormaps.
%
% Example:  To produce a colortable CT of ncol X 3 entries (RGB) of 
%           sequential type and named 'Blues' with 8 colors:
%                   CT=cbrewer('seq', 'Blues', 8);
%           To use this colortable as colormap, simply call:
%                   colormap(CT)
% 
%           To see the various colormaps available according to their types and
%           names, simply call: cbrewer()
%
%  This product includes color specifications and designs developed by
%  Cynthia Brewer (http://colorbrewer.org/).
%
% Author: Charles Robert
% email: tannoudji@hotmail.com
% Date: 06.12.2011
%
% repacked for MS 2.1 - INRA\Olivier Vitrac - 29/01/12 - rev. 13/05/2014

if ~nargin, plot_brewer_cmap(); return, end

% load colorbrewer data
load('colorbrewer.mat')
% initialise the colormap is there are any problems
colormap=[];
if (~exist('interp_method', 'var'))
    interp_method='pchip'; %'cubic';
end

% If no arguments
if (~exist('ctype', 'var') || ~exist('cname', 'var') || ~exist('ncol', 'var'))
    disp(' ')
    disp('INPUT:')
    disp('  - ctype: type of color table *seq* (sequential), *div* (divergent), *qual* (qualitative)')
    disp('  - cname: name of colortable. It changes depending on ctype.')
    disp('  - ncol:  number of color in the table. It changes according to ctype and cname')
    
    disp(' ')
    disp('Sequential tables:')
    z={'Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd',...
             'Purples','RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd'};
    z(:)     
         
    disp('Divergent tables:')
    z={'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn'};
    z(:) 
    
    disp(' ')
    disp('Qualitative tables:')
    %getfield(colorbrewer, 'qual')
    z={'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3'};
    z(:)

    plot_brewer_cmap
    return
end

% Verify that the input is appropriate
ctype_names={'div', 'seq', 'qual'};
if (~ismember(ctype,ctype_names))
    disp('ctype must be either: *div*, *seq* or *qual*')
    colormap=[];
    return
end

if (~isfield(colorbrewer.(ctype),cname))
    disp(['The name of the colortable of type *' ctype '* must be one of the following:'])
    getfield(colorbrewer, ctype)
    colormap=[];
    return
end

if (ncol>length(colorbrewer.(ctype).(cname)))
    disp(' ')
    disp('----------------------------------------------------------------------')
    disp(['The maximum number of colors for table *' cname '* is ' num2str(length(colorbrewer.(ctype).(cname)))])
    disp(['The new colormap will be extrapolated from these ' num2str(length(colorbrewer.(ctype).(cname))) ' values'])
    disp('----------------------------------------------------------------------')
    disp(' ')
    cbrew_init=colorbrewer.(ctype).(cname){length(colorbrewer.(ctype).(cname))};
    colormap=interpolate_cbrewer(cbrew_init, interp_method, ncol);
    colormap=colormap./255;
    return
end

if (isempty(colorbrewer.(ctype).(cname){ncol}))
    
    while(isempty(colorbrewer.(ctype).(cname){ncol}))
        ncol=ncol+1;
    end        
    disp(' ')
    disp('----------------------------------------------------------------------')
    disp(['The minimum number of colors for table *' cname '* is ' num2str(ncol)])
    disp('This minimum value shall be defined as ncol instead')
    disp('----------------------------------------------------------------------')
    disp(' ')
end

colormap=(colorbrewer.(ctype).(cname){ncol})./255;

end

%% external functions
% 
% INTERPOLATE_CBREWER - interpolate a colorbrewer map to ncolors levels
%
% INPUT:
%   - cbrew_init: the initial colormap with format N*3
%   - interp_method: interpolation method, which can be the following:
%               'nearest' - nearest neighbor interpolation
%               'linear'  - bilinear interpolation
%               'spline'  - spline interpolation
%               'cubic'   - bicubic interpolation as long as the data is
%                           uniformly spaced, otherwise the same as 'spline'
%   - ncolors=desired number of colors 
%
% Author: Charles Robert
% email: tannoudji@hotmail.com
% Date: 14.10.2011


function [interp_cmap]=interpolate_cbrewer(cbrew_init, interp_method, ncolors)

% just to make sure, in case someone puts in a decimal
ncolors=round(ncolors);

% How many data points of the colormap available
nmax=size(cbrew_init,1);

% create the associated X axis (using round to get rid of decimals)
a=(ncolors-1)./(nmax-1);
X=round([0 a:a:(ncolors-1)]);
X2=0:ncolors-1;

z=interp1(X,cbrew_init(:,1),X2,interp_method);
z2=interp1(X,cbrew_init(:,2),X2,interp_method);
z3=interp1(X,cbrew_init(:,3),X2, interp_method);
interp_cmap=round([z' z2' z3']);

end

% Plots and identifies the various colorbrewer tables available.
% Is called by cbrewer.m when no arguments are given.
%
% Author: Charles Robert
% email: tannoudji@hotmail.com
% Date: 14.10.2011

function plot_brewer_cmap()
load('colorbrewer.mat')

ctypes={'div', 'seq', 'qual'};
ctypes_title={'Diverging', 'Sequential', 'Qualitative'};
cnames{1,:}={'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn'};
cnames{2,:}={'Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd',...
             'Purples','RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd'};
cnames{3,:}={'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3'};

figure('position', [314 327 807 420])
for itype=1:3
    
    %fh(itype)=figure();
    
    subplot(1,3,itype)
    
    for iname=1:length(cnames{itype,:})
        
        ncol=length(colorbrewer.(ctypes{itype}).(cnames{itype}{iname}));
        fg=1./ncol; % geometrical factor

        X=fg.*[0 0 1 1];
        Y=0.1.*[1 0 0 1]+(2*iname-1)*0.1;
        F=cbrewer(ctypes{itype}, cnames{itype}{iname}, ncol);

        for icol=1:ncol
            X2=X+fg.*(icol-1);
            fill(X2,Y,F(icol, :), 'linestyle', 'none')
            text(-0.1, mean(Y), cnames{itype}{iname}, 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'FontSize', 10, 'FontName' , 'AvantGarde')
            xlim([-0.4, 1])
            hold all
        end % icol
        %set(gca, 'box', 'off')
        title(ctypes_title{itype}, 'FontWeight', 'bold', 'FontSize', 16, 'FontName' , 'AvantGarde')
        axis off
        set(gcf, 'color', [1 1 1])
    end % iname

end %itype

set(gcf, 'MenuBar', 'none')
set(gcf, 'Name', 'ColorBrewer Color maps')
end
