%--------------------------------------------------------------------------------
function [fig,cb] = mapPlot_offline(lon_rho,lat_rho,dat,varargin);
    % -----------------------
    % A way to quickly plot a 2D field
    %
    % Usage:
    % - [fig,cb] = mapPlot(lon,lat,dat,varargin);
    %
    % Inputs:
	% - lon = 2D lon points
	% - lat = 2D lat points
    % - dat = 2D field to plot
    %
    % Varargin:
    % - meta          = option to include structure to name and units (e.g., obj.data.avg.temp)
    % - lonbounds     = x-boundaries (defaults to whole domain)
    % - latbounds     = y-boundaries (defaults to whole domain)
    % - lonticks      = override automatic longitude ticks
    % - latticks      = override automatic latitude ticks
    % - ticks         = 2 = fancy, 1 = on, 0 = off (default set in romsOpt)
    % - background    = background color (default set in romsOpt)
    % - coastcolor    = coast color (default set in romsOpt);
    % - fontsize      = fontsize (default set in romsOpt)
    % - figtype       = 'mfig' (140mm wide), can use 'sfig' (90mm) or 'lfig' (190mm) (default set in romsOpt)
    % - figdim        = figtype multiplier for height (e.g., 1 sets same width and height) 
    % - levels        = hard-coded levels to plot (can't be used with A.prc)
    % - cmap          = colormap(default = thermal for no balance, balance for balance)
    % - caxis         = colorbar limits 
    % - prc           = percentage to limit colorbar axes (if no levels supplied)
    % - bal           = force a balanced colorbar around 0 (1 == yes, 0 == no)
    % - log           = log-scale (1), use with caxis to set limits
    % - XAxisLocation = Override x-axis ticklabels location (default = bottom) 
    % - YAxisLocation = Override y-axis ticklabels location (default = left) 
    % - coast         = m_map map quality (coast, crude, low, high, intermediate)
    % - cblocation:   = location of colorbar (default = eastoutside)  
    % - fontname      = fontname (default = helvetica)
    % -----------------------
    
    % User-inputs
	[nx,ny] = size(dat);
    A.meta          = [];
    A.lonbounds     = [];
    A.latbounds     = [];
    A.lonticks      = [];
    A.latticks      = [];
    A.ticks         = 1;
    A.background    = rgb('LightGray');
    A.coastcolor    = rgb('DimGray');
    A.fontsize      = 6;
    A.figtype       = 'mfig';
    A.coast         = 'coast';
    A.figdim        = round((ny/nx)*10)/10;;
    A.prc           = 2;
    A.bal           = 0;
    A.levels        = [];
    A.cmap          = [];
    A.caxis         = [];
    A.log           = 0;
    A.XAxisLocation = 'bottom';
    A.YAxisLocation = 'left';
	A.cblocation    = 'eastoutside';
    A.fontname      = 'Helvetica';
    A = romsObj.parse_pv_pairs(A,varargin);

    % Make double
    if size(dat)==size(lon_rho);
        lon = double(lon_rho);
        lat = double(lat_rho);
    else
        disp('    ERROR(mapPlot): Check dimensions of input');
		return
    end

    % Get auto-bounds if empty
    if isempty(A.lonbounds) & isempty(A.latbounds);
        lonbounds = [floor(min(lon(:))) ceil(max(lon(:)))];
        latbounds = [floor(min(lat(:))) ceil(max(lat(:)))];
    elseif isempty(A.lonbounds);
        lonbounds = [floor(min(lon(:))) ceil(max(lon(:)))];
        latbounds = A.latbounds;
    elseif isempty(A.latbounds);
        lonbounds = A.lonbounds;
        latbounds = [floor(min(lat(:))) ceil(max(lat(:)))];
    else
        lonbounds = A.lonbounds;
        latbounds = A.latbounds;
    end

    % Set up ticks
    if isempty(A.lonticks)
        dx = round(diff(lonbounds)/60,1)*10;
        lonticks = (lonbounds(1):dx:lonbounds(2));    
        lonticks = round(lonticks);
    else
        lonticks = A.lonticks;
    end
    if isempty(A.latticks)
        dy = round(diff(latbounds)/60,1)*10;
        latticks = (latbounds(1):dy:latbounds(2));
        latticks = round(latticks);
    else
        latticks = A.latticks;
    end

    % Initiate figure
    fig = romsObj.piofigs(A.figtype,A.figdim);

    % Get colormap limits
    if isempty(A.levels)
        clims = romsObj.prclims(dat,'prc',A.prc,'bal',A.bal);
        clevs = linspace(clims(1),clims(2),31); 
    else
        clevs = A.levels;
        clims = [A.levels(1) A.levels(end)];
    end
    if ~isempty(A.caxis);
        clevs(1) = A.caxis(1);
        clevs(end) = A.caxis(2);
        clims = A.caxis;
    end
    if max(dat(:)) == 0 & min(dat(:)) == 0 | isnan(max(dat)) ==1 & isnan(min(dat(:)) == 1);
        dat = nan(size(dat));
        clevs = [0 1];
        clims = linspace(0,1,11);
    else
        dat(dat<clevs(1))   = clevs(1);
        dat(dat>clevs(end)) = clevs(end);
    end

    % Get colormap
    if ischar(A.cmap)
        A.cmap = cmocean(A.cmap,length(clevs)-1);
    elseif isempty(A.cmap)
        if min(clims(:))<0 & max(clims(:)) > 0 | A.bal == 1
            A.cmap = cmocean('balance',length(clevs)-1);
        else
            A.cmap = cmocean('thermal',length(clevs)-1);
        end
    end

    % Make map
    set(0,'CurrentFigure',fig);
    m_proj('mercator','lat',latbounds,'lon',lonbounds); drawnow
    hold on
    if A.ticks == 0
        m_grid('box','on','linestyle','none','xtick',0,'ytick',0,...
               'xticklabels',[],'yticklabels',[],'backgroundcolor',A.background,'fontn',A.fontname); drawnow
    elseif A.ticks == 1
        m_grid('box','on','linestyle','none','xtick',lonticks,...
               'ytick',latticks,'backgroundcolor',A.background,...
			   'fontsize',A.fontsize,'fontn',A.fontname,'yticklabels',latticks,'xticklabels',lonticks);
    elseif A.ticks == 2
        m_grid('box','fancy','linestyle','none','xtick',lonticks,...
               'ytick',latticks,'backgroundcolor',A.background,...
			   'fontsize',A.fontsize,'fontn',A.fontname,'yticklabels',latticks,'xticklabels',lonticks);
    end
    hold on
    m_contourf(lon,lat,dat,clevs,'LineStyle','none');
    cb = colorbar('location',A.cblocation); drawnow
    cb.FontSize = A.fontsize;
    try
        caxis([clims]);
    catch
        caxis([0 1]);
    end
    ax = get(gca);
    if strcmp(A.coast,'coast');
        m_coast('patch',A.coastcolor,'edgecolor','k'); drawnow
    elseif strcmp(A.coast,'crude')
        m_gshhs_c('patch',A.coastcolor,'edgecolor','k'); drawnow
    elseif strcmp(A.coast,'low')
        m_gshhs_l('patch',A.coastcolor,'edgecolor','k'); drawnow
    elseif strcmp(A.coast,'high')
        m_gshhs_h('patch',A.coastcolor,'edgecolor','k'); drawnow
    elseif strcmp(A.coast,'intermediate')
        m_gshhs_i('patch',A.coastcolor,'edgecolor','k'); drawnow
    elseif strcmp(A.coast,'full')
        m_gshhs_f('patch',A.coastcolor,'edgecolor','k'); drawnow
    end
    colormap(gca,A.cmap);
    if A.log == 1 & ~isempty(A.caxis)
        set(gca,'ColorScale','log');
        caxis(A.caxis);
    end
	set(gca,'YAxisLocation',A.YAxisLocation);
	set(gca,'XAxisLocation',A.XAxisLocation);

    % Include meta data if available
    if ~isempty(A.meta)
        title(A.meta.name,'Interpreter','Latex','fontsize',A.fontsize);
        ylabel(cb,A.meta.units,'Interpreter','Latex','fontsize',A.fontsize);
    end

    % Print figure
    if nargout<1
        romsObj.pltjpg(1);
    end

end % end method mapPlot
