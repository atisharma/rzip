% Sets path to use RZIP routines
% cd to the RZIP home directory and execute this .m

rootpath = pwd();

% the RZIP model itself
addpath([rootpath '/rzip'])
addpath([rootpath '/rzip/subroutines'])

% machine description paths
addpath([rootpath '/rzip/tcv'])
addpath([rootpath '/rzip/tcv/data'])
addpath([rootpath '/rzip/tcv-tiles'])

% control synthesis and simulation routines
addpath([rootpath '/control_synthesis'])


% some pre-calculated models & other miscellania
addpath([rootpath '/TCV'])