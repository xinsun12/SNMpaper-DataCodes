function fig = piofigs(pick,ymult)
% Simple function to produce figures of preset width
%
% options:
%	'sfig' = 9cm
%	'mfig' = 14cm
%	'lfig' = 19cm
%   or simply enter a number to set width (e.g., 10 for 10cm)
%
% ymult = number to multiply xwidth by
%

% Return errors
if ymult <= 0
	disp('ymult must be greater than 0');
	return
end

% Grab dim
if ischar(pick)
	if strmatch(pick,'sfig')
		xwidth = 9;
		str    = 'sfig: 9cm';
	elseif strmatch(pick,'mfig');
		xwidth = 14;
		str    = 'mfig: 14cm';
	elseif strmatch(pick,'lfig');
		xwidth = 19;
		str    = 'lfig: 19cm';
	end
else
	xwidth = pick;
end

% Settings (90cm wide)
%disp(['Initializing ',str])
fig             = figure('visible','off');
fig.PaperUnits  = 'Centimeters';
fig.Units       = 'Centimeters';
fig.Position(3) = xwidth;
fig.Position(4) = xwidth*ymult;
fig.Color       = [1 1 1]; 
%disp('Done (invisible)');


