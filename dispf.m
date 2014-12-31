function dispf(varargin)
%DISPF fast wrapper of disp(sprintf(...))
% see help on SPRINTF
% see also FPRINTF (the main difference is that LF is used after disp)

% MS 2.1 - 16/03/08 - INRA\Olivier Vitrac rev. 29/12/12

% Revision history
% 29/12/12 updated help

disp(sprintf(varargin{:})) %#ok<DSPS>