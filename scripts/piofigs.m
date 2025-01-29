function fig = piofigs(pick,ymult)
% pick options:
% sfig = 90mm
% mfig = 140mm
% lfig = 190mm
%
% ymult = number to multiply xwidth by

% Return errors
if ymult <= 0
	disp('ymult must be greater than 0');
	return
end

% Grab dim
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

% Settings (90cm wide)
disp(['Initializing ',str])
fig             = figure('visible','off');
fig.PaperUnits  = 'Centimeters';
fig.Units       = 'Centimeters';
fig.Position(3) = xwidth;
fig.Position(4) = xwidth*ymult;
fig.Color       = [1 1 1]; 
disp('Done (invisible)');


