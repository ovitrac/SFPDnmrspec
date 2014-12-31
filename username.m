function u=username
%USERNAME return the username using the system command whoami

% MS 2.1 - 02/02/08 - INRA\Olivier  - rev. 08/06/09

% revision history
%09/02/08 remove LF...
%08/06/09 use getenv on win boxes

if isunix
    [r,u]=system('whoami');
    u = u((u>='A') & (u<='z'));
else
    u = getenv('USERNAME');
end