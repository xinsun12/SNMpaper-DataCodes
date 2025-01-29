% % % % % % % % % % % % % % %
%  Simple plotting function %
% % % % % % % % % % % % % % %
function [fig,ax,cb] = plotEtaSlice(deg,dep,dat,varargin)
	% Process inputs
	A.levs          = [linspace(-0.05,0.05,52)];
	A.cmap          = [cmocean('balance',51)]; 
	A.YAxisLocation = 'left';
	A.hatch         = 0;
	A.hatchthresh   = 0;
	A.hatchdata     = [];
	A.YLim          = [-1000 0];
	A.XLim          = [min(deg(:)) max(deg(:))];
	A.cblocation    = 'eastoutside';
	A.figtype       = 'mfig';
	A.figdim        = 0.33;
	A.fontsize      = [];
	A.XTick         = [];
	A = parse_pv_pairs(A,varargin);

	% Get fig
	fig = piofigs(A.figtype,A.figdim);

	% Plot contour
	dat(dat<A.levs(1)) = A.levs(1);
	dat(dat>A.levs(end)) = A.levs(end);
	[c1,h1] = contourf(deg,dep,dat,A.levs,'linestyle','none');
	hold on
	ax = gca;
	xlim([A.XLim]);	
	ylim([A.YLim]);
	xlabel('Longitude','Interpreter','Latex');
	if isempty(A.XTick);
		xtck = get(ax,'XTick');
		xlbl = get(ax,'XTickLabel');
	else
		xtck = A.XTick;
		set(ax,'XTick',xtck);
		set(ax,'XTickLabel',{xtck});
		xlbl = get(ax,'XTickLabel');
	end
	for i = 1:length(xlbl)
		if str2num(xlbl{i})>180
			%newlbl{i} = [num2str(abs(str2num(xlbl{i})-360)),char(176),'W'];
			newlbl{i} = [num2str(abs(str2num(xlbl{i})-360)),'$^o$','W'];
		else
			newlbl{i} = [num2str(str2num(xlbl{i})),'$^o$','E'];
		end
	end
	set(ax,'XTick',xtck);
	set(ax,'XTickLabel',newlbl);
	ylabel('Depth (m)','Interpreter','Latex');
	cb = colorbar('location',A.cblocation);
	caxis([A.levs(1) A.levs(end)]);
	set(ax,'ColorMap',A.cmap);
	set(ax,'Color',rgb('DimGray'));
	set(ax,'YAxisLocation',A.YAxisLocation)
	set(ax,'TickLabelInterpreter','Latex');
	set(fig,'inverthardcopy','off');
	if ~isempty(A.fontsize)
		set(ax,'FontSize',A.fontsize);
	end
end
