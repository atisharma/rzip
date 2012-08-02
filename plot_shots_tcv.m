% plot figures for various TCV shots

Machine_name = 'tcv';

time = .3;
shots = [13333 16685 5650 6479 9490 10396 13328 13331 13332 13334 14444];

for shot = shots
	plot_machine(Machine_name, 1, shot, time);
	axis([0 2.5 -2 1.5]);
	eval(['print -dpng ./temp/tcv_' num2str(shot) ]);
end