function f = PlotMuscleGeometry(R, info)
%Plot muscle geometry 

% Info on figure
lw = 1.5; % linewidth
f  = figure('Name','Muscle Geometry');

color_ext = [148 93 94]./255;
color_flex= [221 167 123]./255;

nPhases = 3; 
if nPhases > 1
    fTabGroup = uitabgroup;
end

% Cs = linspecer(3);

for i=1:nPhases
    % set the name of the tab
    if nPhases>1
        file = {'Length','Force','LMT-MA'};
        tab = uitab(fTabGroup, 'Title', file{i});
        axes('parent',tab);
    end

    if i ==1
        subplot(221)
        plot(R.C.lM_ext,'Color',color_ext,'LineWidth',lw); hold on; 
        plot(R.C.lM_flex,'Color',color_flex','LineWidth',lw); hold on
        title('lM'); ylabel('[m]'); box off; 
        subplot(222)
        plot(R.lMtilda(1,:),'Color',color_ext,'LineWidth',lw); hold on
        plot(R.lMtilda(2,:),'Color',color_flex,'LineWidth',lw); hold on
        title('lMtilda'); box off
        subplot(223)
        plot(R.C.lT_ext,'Color',color_ext,'LineWidth',lw); hold on;
        plot(R.C.lT_flex,'Color',color_flex,'LineWidth',lw); hold on;
        title('lT'); ylabel('[m]'); xlabel('frames (200Hz)'); box off; legend({'Ext','Flex'},'Location','Best'); 
        subplot(224)
        plot(R.C.lTtilda_ext,'Color',color_ext,'LineWidth',lw); hold on;
        plot(R.C.lTtilda_flex,'Color',color_flex,'LineWidth',lw); hold on;
        title('lTtilda'); xlabel('Frames (200Hz)'); box off; 
        %suptitle([info.subj ' Trial ', num2str(info.trial)])
        
    elseif i == 2
        subplot(231)
        plot(R.C.Fce_ext,'Color',color_ext,'LineWidth',lw); hold on
        plot(R.C.Fce_flex,'Color',color_flex,'LineWidth',lw); hold on
        title('Fce'); box off; 
        subplot(232)
        plot(R.C.Fpe_ext,'Color',color_ext,'LineWidth',lw); hold on
        plot(R.C.Fpe_flex,'Color',color_flex,'LineWidth',lw); hold on
        title('Fpe'); box off; legend({'Ext','Flex'},'Location','Best'); 
        subplot(233)
        plot(R.C.fse_ext,'Color',color_ext,'LineWidth',lw); hold on
        plot(R.C.fse_flex,'Color',color_flex,'LineWidth',lw); hold on
        title('Fse'); box off; 
        subplot(234)
        plot(R.C.FM_ext,'Color',color_ext,'LineWidth',lw); hold on
        plot(R.C.FM_flex,'Color',color_flex,'LineWidth',lw); hold on
        title('FM'); xlabel('Frames (200Hz)'); box off;
        subplot(235)
        plot(R.C.FT_ext,'Color',color_ext,'LineWidth',lw); hold on
        plot(R.C.FT_flex,'Color',color_flex,'LineWidth',lw); hold on
        title('FT'); xlabel('Frames (200Hz)'); box off
        %suptitle([info.subj ' Trial ', num2str(info.trial)])
        
    else
        subplot(211)
        plot(R.C.lMT_ext,'Color',color_ext,'LineWidth',lw); hold on
        plot(R.C.lMT_flex,'Color',color_flex,'LineWidth',lw); hold on
        title('lMT'); box off;
        subplot(212)
        plot(R.C.MA_ext,'Color',color_ext,'LineWidth',lw); hold on
        plot(R.C.MA_flex,'Color',color_flex,'LineWidth',lw); hold on
        title('MA'); box off; legend({'Ext','Flex'},'Location','Best'); 
    end
        
end

