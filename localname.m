function name = localname
% LOCALNAME return the netbios name of the local machine

% MS 1.0 - 26/10/04 - INRA\Olivier Vitrac - rev. 03/02/08
%
% revision history
% 03/02/08 unix compatibility

if isunix
    [r,name] = system('hostname');
else
    name = getenv('COMPUTERNAME');
end
name = name(name>32);