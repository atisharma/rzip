%
%    RZIP - a linear tokamak plasma equilibrium response model
%    
%    Copyright (C) 2000 J Wainwright, A Sharma
%    ati.sharma@gmail.com
%    
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
%
%    In addition, if you use the code or a derivative, please reference the
%    appropriate papers.
%
% This demo illustrates some features of the RZIP model and its routines
% using TCV shot 13333 as an example
% (the shot used for Coutlis' ID work and the control work)
%
% A Sharma 2000

clear all

% must be in this directory to set the paths
set_rzip_paths; cd temp;

% Shot details
Machine_name = 'tcv';
shot = 13333;
time = .3;

% The 18 coils of TCV - the 19th `fast' G-coil is left open-circuit for this shot
Used_coils = 1:18;

% The number of modes in the passive structure reduction
N_eigen_modes_array = 38; 

% [Plasma resistance, d(Plasma resistance)/dr, d(beta_p)/dt, d(l_i)/dt] - these are default values and can be 'fitted' to the data in a grey-box way - see NF JT-60U paper
plasma_inputs = [0 0 0 0]; 

%---------------------------------------------------------------

% plot the machine description
plot_machine(Machine_name, 4, shot, time)

% Make the RZIP model (version 3)
[A, B, C, D, curlyM, curlyR, plasma_outputs, IFLAG] = rzip_v3(shot, time, Used_coils, N_eigen_modes_array, plasma_inputs, Machine_name, 0, '');

% remove zeros eigenvalues
[A,B,C,D] = schurred(A,B,C,D);

% The control parameters
tcv_observer
C = Cmod(1:5,:)*C(1:94,:);
D = Dmod(1:5,:)*D(1:94,:);

% for ss object
G = ss(A,B,C,D);

% create and save power supply model
pson = 1; coils=1:18;
tcv_powersupplies

% Include p/s in model
G=G(:,1:18)*Gps;

% synthesise reduced-order NLCF H-inf controller using weights
hinf_control_synthesis

% prepare a simulation
makesim

% simulate with H-inf controller
runsim

% load PID controller for shot 13333 & get in useful form
PID_c2d; PIDforsim;

% simulate with PID controller & compare
runsim_PID

cd ..

% print out
figure(1);
print -dpng ./temp/RZIP_controlled_sim




