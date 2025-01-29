function [patch] = shaded_1d(lowX,highX,Y,ec,fc,fa);
% Function to make shaded 1D patches
% Usage:
% - [patch] = shaded_1d(lowX,highX,Y,meanX,ec,fc,fa);
%
% Inputs:
% - lowX  = lower bound (min? mean - 2std? etc)
% - highX = upper bound (max? mean + 2std? etc) 
% - Y     = Y values (depth? etc) 
%
% Optional Inputs:
% - ec = 'edgecolor', default = none
% - fc = 'facecolor', default = rgb('DimGray')
% = fa = 'facealpha', default = 0.7
if nargin < 3
	disp('error: see shaded_1d, need 3 inputs');
	return
elseif nargin < 4
	ec = 'none';
	fc = rgb('DimGray');
	fa = 0.7;
end
x1  = [lowX];
y1  = [Y];
dat = [x1+y1];
x1  = [x1(~isnan(dat))];
y1  = [y1(~isnan(dat))];
% Right side
x2  = [highX];
y2  = [Y];
dat = [x2+y2];
x2  = [x2(~isnan(dat))];
y2  = [y2(~isnan(dat))];
x2  = flipud(x2);
y2  = flipud(y2);
% Connect
X = [x1 ;x2];
Y = [y1 ;y2];
patch = fill(X,Y,fc);
patch.EdgeColor = ec;
patch.FaceAlpha = fa;
