% Script to compare output from 2 runs
% One with nitrite-oxidizers (obj1)
% One with comammox (obj2)

% Switch needed to process output
if (0)
	clearvars -except obj1 obj2
	warning off
	addpath /data/project1/demccoy/romsObj/
	addpath /data/project1/demccoy/ROMS/peru_chile_0p1/analysis/microbes/paper

	% Make plot directory
	mkdir plots

	% List runs
	simName  = 'peru_chile_0p1';
	runName1 = 'microbes_eth_obligate_baseline_Km';
	runName2 = 'microbes_eth_obligate_baseline_Km_comammox';
	rTit1    = 'obligate (nitrox)';
	rTit2    = 'obligate (comammox)';
	domain   = [41 401 261 461];

	% Get objects
	try; obj1;
	catch
		obj1 = initROMS(romsObj,simName,runName1,'domain',domain);
	end
	try; obj2;
	catch
		obj2 = initROMS(romsObj,simName,runName2,'domain',domain);
	end
else
	addpath scripts
end

% Get plotdir
pltdir = ['plots/'];

% Extract profiles from OMZ
if (0)
	% Options
	vars = {'O2','DOC','DOCR','NO2',...
	        'AER','NAR','NAI','NAO',...
			'NIR','NIO','NOS',...
			'AOA','NOB','AOX',...
			'DENITRIF1','DENITRIF2','DENITRIF3',...
			'DENITRIF4','DENITRIF5','DENITRIF6',...
			'AMMOX','NITROX','ANAMMOX',...
			'SP_NO2_UPTAKE','DIAT_NO2_UPTAKE','DIAZ_NO2_UPTAKE'};
	prof_lon = [226:231];
	prof_lat = [117:122];
	zgrid = [-805:1:0];
	files1 = 1:240;
	files2 = 1:240;

	% Initialize output matrices
	for i = 1:length(vars)
		OUT1.(vars{i}).mean = nan(length(zgrid),length(files1)+1);
		OUT1.(vars{i}).std  = nan(length(zgrid),length(files1)+1);
		OUT2.(vars{i}).mean = nan(length(zgrid),length(files1)+1);
		OUT2.(vars{i}).std  = nan(length(zgrid),length(files1)+1);
	end
		
	% Cycle through variables
    for i = 1:length(vars)

		% Write depths
		OUT1.(vars{i}).mean(:,1) = flipud(-zgrid');
		OUT2.(vars{i}).mean(:,1) = flipud(-zgrid');
		OUT1.(vars{i}).std(:,1)  = flipud(-zgrid');
		OUT2.(vars{i}).std(:,1)  = flipud(-zgrid');

		% Interpolate to zgrid
		for t = 1:length(files1)
			% Initialize profile counter
			profcnt = 1;

			% Load depth and data
			obj1 = clearROMS(obj1); 
			obj1 = loadData(obj1,vars(i),files1(t));
			obj1 = loadDepth(obj1,files1(t));
			obj2 = clearROMS(obj2); 
			obj2 = loadData(obj2,vars(i),files1(t));
			obj2 = loadDepth(obj2,files1(t));
			% Interpolate to zgrid
			for j = 1:length(prof_lon)
				for k = 1:length(prof_lat)
					tmp.data = squeeze(obj1.data.avg.(vars{i}).data(prof_lon(j),prof_lat(k),:));
					tmp.z_r  = squeeze(obj1.grid.avg.z_r(prof_lon(j),prof_lat(k),:));
					nanidx   = tmp.data + tmp.z_r;
					if i > 14
						tmp.data = tmp.data .* 86400;
						if ismember(vars(i),{'DENITRIF3','ANAMMOX'})
							tmp.data = tmp.data .* 2;
						end
					end
					out1.(vars{i})(:,profcnt) = interp1(tmp.z_r(~isnan(nanidx)),tmp.data(~isnan(nanidx)),zgrid);
					tmp.data = squeeze(obj2.data.avg.(vars{i}).data(prof_lon(j),prof_lat(k),:));
					tmp.z_r  = squeeze(obj2.grid.avg.z_r(prof_lon(j),prof_lat(k),:));
					nanidx   = tmp.data + tmp.z_r;
					if i > 14
						tmp.data = tmp.data .* 86400;
						if ismember(vars(i),{'DENITRIF3','ANAMMOX'})
							tmp.data = tmp.data .* 2;
						end
					end
					out2.(vars{i})(:,profcnt) = interp1(tmp.z_r(~isnan(nanidx)),tmp.data(~isnan(nanidx)),zgrid);
					profcnt = profcnt + 1;
				end
			end
			OUT1.(vars{i}).mean(:,t+1) = flipud(nanmean(out1.(vars{i}),2)); 
			OUT2.(vars{i}).mean(:,t+1) = flipud(nanmean(out2.(vars{i}),2));
			OUT1.(vars{i}).std(:,t+1)  = flipud(nanstd(out1.(vars{i}),0,2)); 
			OUT2.(vars{i}).std(:,t+1)  = flipud(nanstd(out2.(vars{i}),0,2));
		end
	end

	% Fill missing values with nearest neighbor
	for i = 1:length(vars)
		for j = 1:size(OUT1.(vars{i}).mean,2);
			OUT1.(vars{i}).mean(:,j) = fillmissing(OUT1.(vars{i}).mean(:,j),'nearest');
			OUT2.(vars{i}).mean(:,j) = fillmissing(OUT2.(vars{i}).mean(:,j),'nearest');
			OUT1.(vars{i}).std(:,j)  = fillmissing(OUT1.(vars{i}).std(:,j),'nearest');
			OUT2.(vars{i}).std(:,j)  = fillmissing(OUT2.(vars{i}).std(:,j),'nearest');
		end
	end

	% Write to matrix
	for i = 1:length(vars)
		cmd = ['rm nitrox.xls']; system(cmd);
		cmd = ['rm comammox.xls']; system(cmd);
		writematrix(OUT1.(vars{i}).mean,'nitrox.xls','FileType','spreadsheet','Sheet',[vars{i},'_mean']);
		writematrix(OUT2.(vars{i}).mean,'comammox.xls','FileType','spreadsheet','Sheet',[vars{i},'_mean']);
	end
end

% Extract and plot transects of O2, NO2, and NO2-->N2 denitrfier biomass from OMZ
if (1)
	% Options
	vars = {'O2','NO2','NIO'};
	eta_ind = [120];
	zgrid = [-805:1:0];
	files1 = 181:240;
	files2 = 181:240;
	dep_lim = [0 600];
	fontsize = 6;
    prof_lon = -83+360;

	% Process or just load data directly
	if (0)
		% Slice output
		obj1 = sliceROMS(obj1,vars,'eta',eta_ind,files1);
		obj2 = sliceROMS(obj2,vars,'eta',eta_ind,files2);

		% Loop through vars, make slice, plot
		for i = 1:length(vars)
			% Call zslice
			tmpdat1 = squeeze(obj1.data.avg.(vars{i}).slice);
			tmpdat1(tmpdat1<0) = 0;
			tmpdat2 = squeeze(obj2.data.avg.(vars{i}).slice);
			tmpdat2(tmpdat2<0) = 0;
			[a,b,c] = size(tmpdat1);
			clear tmp_dat1 tmp_dat2
			% Interpolate to zdeps
			for j = 1:a
				for k = 1:c
					tmpdep1 = squeeze(obj1.slice.depth(j,:,1,k));
					tmpdep2 = squeeze(obj2.slice.depth(j,:,1,k));
					ind1 = find(~isnan(tmpdat1(j,:,k))==1);
					ind2 = find(~isnan(tmpdat2(j,:,k))==1);
					if ~isempty(ind1)
						tmp_dat1(j,:,k) = interp1(tmpdep1,tmpdat1(j,:,k),zgrid);
					else
						tmp_dat1(j,:,k) = nan(size(zgrid));
					end
					if ~isempty(ind2)
						tmp_dat2(j,:,k) = interp1(tmpdep2,tmpdat2(j,:,k),zgrid);
					else
						tmp_dat2(j,:,k) = nan(size(zgrid));
					end
				end
			end
			% Save averages
			out1.(vars{i}) = nanmean(tmp_dat1,3);
			out2.(vars{i}) = nanmean(tmp_dat2,3);
		end

		% Get grid
		tmpdeg = obj1.slice.deg(:,1,1,1);
		tmpdeg = repmat(tmpdeg,[1 length(zgrid)]);
		tmpdep = zgrid';
		tmpdep = repmat(tmpdep,[1 size(tmpdeg,1)])';

		% Save output for Xin?
		save('data/Sun_et_al_model_data1.mat','out1','out2','tmpdeg','tmpdep');
	else
		load('data/Sun_et_al_model_data1.mat');
	end

	% Get data
	dat = out1.NO2 - out2.NO2;
	lvls = linspace(-3,3,22);
	cmap = cmocean('balance',21);
	cmap(11,:) = [1 1 1];

	% Make figure
	[fig,ax,cb] = plotEtaSlice(tmpdeg,-tmpdep,dat,...
		'levs',lvls,...
		'cmap',cmap,...
		'YLim',[dep_lim],...
		'XLim',[256.7489 282],...
		'XTick',[260 270 280],...
		'figtype','sfig',...
		'figdim',0.33,...
		'cblocation','eastoutside',...
		'fontsize',fontsize);
        hold on
	contour(tmpdeg,-tmpdep,out1.O2,[1 1],...
		'linestyle','--','color',rgb('DimGray'),'linewidth',0.25);
	ylim([dep_lim]);
	plot(prof_lon,dep_lim(end),'ro','markersize',2,'markerfacecolor','r');
	ylabel(cb,'$\mu$M','Interpreter','Latex');
	cb.TickLabelInterpreter = 'Latex';
	load('data/transect_position_east.mat');
	ax.Position = axpos;
	cb.Position = cbpos;
	ax.YDir = 'Reverse';
	cb.Ticks = [-3 -2 -1 0 1 2 3];
	export_fig('-png',['plots/no2_trans_cmp'],'-m5');
	close all

	% Get data
	dat = out1.NIO - out2.NIO;
	lvls = linspace(-0.2,0.2,22);
	cmap = cmocean('balance',21);
	cmap(11,:) = [1 1 1];

	% Make figure
	[fig,ax,cb] = plotEtaSlice(tmpdeg,-tmpdep,dat,...
		'levs',lvls,...
		'cmap',cmap,...
		'YLim',[dep_lim],...
		'XLim',[256.7489 282],...
		'XTick',[260 270 280],...
		'figtype','sfig',...
		'figdim',0.33,...
		'cblocation','eastoutside',...
		'fontsize',fontsize);
        hold on
	contour(tmpdeg,-tmpdep,out1.O2,[1 1],...
		'linestyle','--','color',rgb('DimGray'),'linewidth',0.25);
	ylim([dep_lim]);
	ylabel(cb,'$\mu$M','Interpreter','Latex');
	cb.TickLabelInterpreter = 'Latex';
	load('data/transect_position_east.mat');
	ax.Position = axpos;
	cb.Position = cbpos;
	ax.YDir = 'Reverse';
	export_fig('-png',['plots/nio_trans_cmp'],'-m5');
	close all

    % Get data
    dat = out1.NO2;
    lvls = linspace(0,3,21);
    cmap = cmocean('tempo',20);
    cmap(1,:) = [1 1 1];

    % Make figure
    [fig,ax,cb] = plotEtaSlice(tmpdeg,-tmpdep,dat,...
        'levs',lvls,...
        'cmap',cmap,...
        'YLim',[dep_lim],...
        'XLim',[256.7489 282],...
        'XTick',[260 270 280],...
        'figtype','sfig',...
        'figdim',0.33,...
        'cblocation','eastoutside',...
        'fontsize',fontsize);
        hold on
    contour(tmpdeg,-tmpdep,out1.O2,[1 1],...
        'linestyle','--','color',rgb('DimGray'),'linewidth',0.25);
    ylim([dep_lim]);
	plot(prof_lon,dep_lim(end),'ro','markersize',2,'markerfacecolor','r');
    ylabel(cb,'$\mu$M','Interpreter','Latex');
    cb.TickLabelInterpreter = 'Latex';
    load('data/transect_position_east.mat');
    ax.Position = axpos;
    cb.Position = cbpos;
    ax.YDir = 'Reverse';
    export_fig('-png',['plots/no2_trans_nitrox'],'-m5');
    close all
	
end

% Extract and plot NO2 and production/consumption profiles from OMZ
if (1)

	% Options
	vars = {'O2','NO2','DENITRIF1','DENITRIF2','DENITRIF5','NITROX','NO2AMMOX','ANAMMOX','NO2_UPTAKE','O2_PRODUCTION','O2_CONSUMPTION'};
	fact = [  1     1       1           -1          -1        -1        1          -1         -1           1                -1        ];
    prof_lon = [226:231];
    prof_lat = [117:122];
    zgrid = [-805:1:0];
	files1 = 181:240;
	fontsize = 6;

	% List all rates from model
	rates = {'OXIC_REMIN','DENITRIF1','DENITRIF4','DENITRIF6','DENITRIF2','DENITRIF5','DENITRIF3',...
             'AMMOX','NITROX','ANAMMOX','NO2AMMOX','N2OAMMOX','NO3_UPTAKE','NO2_UPTAKE'};
	rtits = {'$R^{aer}_{oxic}','$R^{nar}_{den1}$','$R^{nai}_{den4}$',...
             '$R$^{nao}_{den6}$','$R^{nir}_{den2}$','$R^{nio}_{den5}$','$R^{nos}_{den6}$',...
             '$R^{aoa}_{ao}$','$R^{nob}_{no}$','$R^{aox}_{ax}$','$R^{aoa}_{ao,no2}$','$R^{aoa}_{ao,n2o}$',...
             '$R_{V,no3}$','$R_{V,no2}$'};
	rclrs = colormix(length(rates)-4,{'k','w','y'});
	rclrs(end+1,:) = rclrs(8,:); % copy color for NO2AMMOX
	rclrs(end+1,:) = rclrs(8,:); % copy color for N2OAMMOX
	rclrs(end+1,:) = rgb('Black'); % Uptake to black
	rclrs(end+1,:) = rgb('Black'); % Uptake to black

	% Process data or load directly?
	if (0)
		% Grab profiles	
		obj1 = loadDepth(obj1,files1);
		for i = 1:length(vars)
			profcnt = 1;
			% Load data
			obj1 = loadData(obj1,vars(i),files1);
			% Interpolate to zgrid
			for j = 1:length(prof_lon)
				for k = 1:length(prof_lat)
					for t = 1:12
						tmp.data = squeeze(obj1.data.avg.(vars{i}).data(prof_lon(j),prof_lat(k),:,t)).*fact(i);
						if i > 2
							tmp.data = tmp.data .* 86400; % convert to per day
						end
						tmp.z_r  = squeeze(obj1.grid.avg.z_r(prof_lon(j),prof_lat(k),:,t));
						nanidx   = tmp.data + tmp.z_r;
						tmp.out.(vars{i})(:,profcnt) = interp1(tmp.z_r(~isnan(nanidx)),tmp.data(~isnan(nanidx)),zgrid);
						profcnt = profcnt + 1;
					end
				end
			end
		end

		% Compute SMS
		tmp.out.NO2_SMS = tmp.out.DENITRIF1 + tmp.out.DENITRIF2 + tmp.out.DENITRIF5  + tmp.out.NITROX + ...
						  tmp.out.NO2AMMOX  + tmp.out.ANAMMOX   + tmp.out.NO2_UPTAKE; 
		tmp.out.O2_SMS  = tmp.out.O2_PRODUCTION + tmp.out.O2_CONSUMPTION;
		vars{end+1} = 'NO2_SMS';
		vars{end+1} = 'O2_SMS';

		% Compute statistics
		for i = 1:length(vars)
			if (0)
				[percs] = prctile(tmp.out.(vars{i}),[25 50 75],2);
				out.(vars{i}).p25 = percs(:,1);
				out.(vars{i}).p50 = percs(:,2);
				out.(vars{i}).p75 = percs(:,3);
			else
				tmpmean = nanmean(tmp.out.(vars{i}),2);
				tmpstd = nanstd(tmp.out.(vars{i}),0,2);
				out.(vars{i}).p25 = tmpmean - 2*tmpstd;
				out.(vars{i}).p50 = tmpmean;
				out.(vars{i}).p75 = tmpmean + 2*tmpstd;
			end
		end

		save('data/Sun_et_al_model_data2.mat','out','zgrid');
	else
		load('data/Sun_et_al_model_data2.mat');
		vars{end+1} = 'NO2_SMS';
		vars{end+1} = 'O2_SMS';
	end

	% Tracer and O2 subplot
	% Initialize figure 
	fig = piofigs('sfig',1);
	[ax,h1,h2] = plotxx(nan(size(zgrid)),zgrid,nan(size(zgrid)),zgrid,...
		{['O$_2$ ($\mu$M)'],['NO$^{-}_2$ ($\mu$M)']},...
		{'Depth (m)','Depth (m)'});
	% Plot O2
	set(fig,'CurrentAxes',ax(1));
	hold on
	[patch1] = shaded_1d(out.O2.p25,out.O2.p75,-zgrid','none',rgb('DimGray'),0.5);
	plot(out.O2.p50,-zgrid,'-','color',rgb('Black'),'linewidth',1);
	ax(1).XLabel.Interpreter = 'Latex';
	ax(1).YLabel.Interpreter = 'Latex';
	ax(1).XLim = [0 250]; 
	ax(1).YLim = [0 600];
	% Plot NO3
	set(fig,'CurrentAxes',ax(2));
	hold on
	[patch2] = shaded_1d(out.NO2.p25,out.NO2.p75,-zgrid','none',rgb('Green'),0.5);
	plot(out.NO2.p50,-zgrid,'-','color',rgb('DarkGreen'),'linewidth',1);
	ax(2).XLabel.Interpreter = 'Latex';
	ax(2).YLabel.Interpreter = 'Latex';
	ax(2).XLim = [0 2]; 
	ax(2).YLim = [0 600];
    ax(2).YTickLabel = [];
    ax(2).YLabel.String = [];
    ax(2).XColor = rgb('DarkGreen');
	grid on
	ax(1).YDir = 'reverse';
	ax(2).YDir = 'reverse';
    ax(1).Position(4) = ax(1).Position(4)*0.8;
    ax(2).Position(4) = ax(2).Position(4)*0.8;
    ax(1).Position(3) = ax(1).Position(3)*0.6;
    ax(2).Position(3) = ax(2).Position(3)*0.6;
    ax(1).Position(2) = ax(1).Position(2)+0.05;
    ax(2).Position(2) = ax(2).Position(2)+0.05;
	ax(1).Position(1) = ax(1).Position(1)+0.10;
	ax(2).Position(1) = ax(2).Position(1)+0.10;
	ax(1).FontSize = fontsize; 
	ax(2).FontSize = fontsize; 
    axpos = ax(1).Position;
	% Save
	export_fig('-png',[pltdir,'NO2_vs_O2_prof'],'-m5');
	close all

	% O2 Sources & sinks subplot
	fig = piofigs('sfig',1);
	ax = gca;
	hold on
	lcnt = 1;
	vars = {'O2_PRODUCTION','O2_CONSUMPTION','O2_SMS'};
	clrs = [rgb('FireBrick');rgb('RoyalBlue');rgb('Black')];
	% Plot dummy legend
	plot(NaN,NaN,'-','color',clrs(1,:));
	plot(NaN,NaN,'-','color',clrs(2,:));
	plot(NaN,NaN,'-','color',clrs(3,:));
	legend_entries{1} = '$R^{prod}_{o2}$';
	legend_entries{2} = '$R^{cons}_{o2}$';
	legend_entries{3} = '$J_{o2}$';
	[l] = legend(legend_entries,'location','southeast','Interpreter','Latex');
	l.Box = 'off';
	l.AutoUpdate = 'off';
	l.ItemTokenSize = [10 6];	
	for i = 1:length(vars)
		hold on
		[patch2(i)] = shaded_1d(out.(vars{i}).p25,out.(vars{i}).p75,-zgrid','none',clrs(i,:),0.5);
		hold on
		plot(out.(vars{i}).p50,-zgrid,'-','color',clrs(i,:),'linewidth',1);
	end
	xlabel('mmol O$_2$ m$^{-3}$ d$^{-1}$','Interpreter','Latex');
	set(ax,'YTickLabel',[]);
	title('O$_2$: Sources and Sinks','Interpreter','Latex');
	grid on
	ax.Position = axpos;
	set(ax,'YDir','Reverse');
	xlim([-3 3]);
	ylim([0 600]);
	ax.Position = axpos;
	ax.FontSize = fontsize; 
	% Save
    export_fig('-png',[pltdir,'O2_sms'],'-m5');
    close all

	% NO2 Sources & sinks subplot
	fig = piofigs('sfig',1);
	ax = gca;
	hold on
	lcnt = 1;
	vars = {'DENITRIF1','DENITRIF2','DENITRIF5','NITROX','NO2AMMOX','ANAMMOX','NO2_UPTAKE','NO2_SMS'};
	% Plot dummy legend
	lcnt = 1;
	for i = 1:length(vars)
		idx = find(strcmp(vars{i},rates)==1);
		if ~isempty(idx)
			plot(NaN,NaN,'-','color',rclrs(idx,:));
			legend_entries{lcnt} = rtits{idx};
		else
			plot(NaN,NaN,'-','color',rgb('Black'));
			legend_entries{lcnt} = '$J_{no2}$';
		end
		lcnt = lcnt + 1;
	end
	[l] = legend(legend_entries,'location','southeast','Interpreter','Latex');
	l.Box = 'off';
	l.AutoUpdate = 'off';
	l.ItemTokenSize = [10 6];	
	for i = 1:length(vars)
		idx = find(strcmp(vars{i},rates)==1);
		if ~isempty(idx)
			[patch2(i)] = shaded_1d(out.(vars{i}).p25,out.(vars{i}).p75,-zgrid','none',rclrs(idx,:),0.5);
			plot(out.(vars{i}).p50,-zgrid,'-','color',rclrs(idx,:),'linewidth',1);
		else
			[patch2(i)] = shaded_1d(out.(vars{i}).p25,out.(vars{i}).p75,-zgrid','none',rgb('Black'),0.5);
			plot(out.(vars{i}).p50,-zgrid,'-','color',rgb('Black'),'linewidth',1);
		end
	end
	xlabel('mmol N m$^{-3}$ d$^{-1}$','Interpreter','Latex');
	set(ax,'YTickLabel',[]);
	title('NO$^{-}_2$: Sources and Sinks','Interpreter','Latex');
	grid on
	ax.Position = axpos;
	set(ax,'YDir','Reverse');
	xlim([-0.06 0.06]);
	ylim([0 600]);
	ax.FontSize = fontsize; 
	% Save
    export_fig('-png',[pltdir,'NO2_sms'],'-m5');
    close all

	% Fraction of production vs consumption
	prod_tot = out.DENITRIF1.p50 + out.NO2AMMOX.p50;
	cons_tot = out.NITROX.p50 + out.DENITRIF2.p50 + out.DENITRIF5.p50 + out.ANAMMOX.p50 + out.NO2_UPTAKE.p50;
	pvars = {'DENITRIF1','NO2AMMOX'};
	for i = 1:length(pvars)
		prod_frac(:,i) = 100*out.(pvars{i}).p50./prod_tot;
		idx = find(strcmp(pvars{i},rates)==1);
		pclrs(i,:) = rclrs(idx,:);
	end
	cvars = {'DENITRIF2','DENITRIF5','NITROX','ANAMMOX','NO2_UPTAKE'};
	for i = 1:length(cvars)
		cons_frac(:,i) = -100*out.(cvars{i}).p50./cons_tot;
		idx = find(strcmp(cvars{i},rates)==1);
		cclrs(i,:) = rclrs(idx,:);
	end
	fig = piofigs('sfig',1);
	ax = gca;
	hold on
	pb = barh(-zgrid,prod_frac,'stacked');
	cb = barh(-zgrid,cons_frac,'stacked');
	for i = 1:length(pb)
		pb(i).EdgeColor = 'none';
		pb(i).FaceColor = pclrs(i,:);
		pb(i).BarWidth = 1;
	end
	for i = 1:length(cb)
		cb(i).EdgeColor = 'none';
		cb(i).FaceColor = cclrs(i,:);
		cb(i).BarWidth = 1;
	end
	xlim([-100 100]);
	ylim([0 600]);
	set(ax,'YTickLabel',[]);
	set(ax,'YDir','Reverse');
	ax.Position = axpos;
	grid on
	title('NO$^{-}_2$: Consumption vs Production','Interpreter','Latex');
	xlabel('Fraction ($\%$)','Interpreter','Latex');
	ax.FontSize = fontsize; 
	% Save
    export_fig('-png',[pltdir,'NO2_sms_frac'],'-m5');
    close all

	try	
		obj1 = clearROMS(obj1);
	catch
	end
end
