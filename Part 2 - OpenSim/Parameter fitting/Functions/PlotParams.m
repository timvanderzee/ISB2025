function g = PlotParams(R,info)
% Plot params of tracking simulations

% Info on figure
g     = figure('Name','Optimized Parameters');
color = [0.6 0.6 0.6]; 
color_ext = [148 93 94]./255;
color_flex= [221 167 123]./255;

subplot(141)
bar(1,R.a(1),'FaceColor',color_ext,'EdgeColor',color_ext); hold on
bar(2,R.a(2),'FaceColor',color_flex,'EdgeColor',color_flex); hold on
line([0 3],[R.bounds.Ub.a R.bounds.Ub.a],'Color','k','LineWidth',1.5);hold on %Ub
line([0 3],[R.bounds.Lb.a R.bounds.Lb.a],'Color','k','LineWidth',1.5);hold on %Lb
title('a'); xticks([1 2]); xticklabels({'Ext','Flex'})
%title([info.subj ' Trial ', num2str(info.trial)])
box off; 

subplot(142)
bar(1,R.kFpe,'FaceColor',color,'EdgeColor',color); hold on
line([0 2],[R.bounds.Ub.kFpe R.bounds.Ub.kFpe],'Color','k','LineWidth',1.5); hold on; %Ub
line([0 2],[R.bounds.Lb.kFpe R.bounds.Lb.kFpe],'Color','k','LineWidth',1.5); hold on; %Lb
line([0 2],[0.1 0.1],'Color','r','LineWidth',1.5); hold on; % nominal
ylim([0 0.2]);
title('kFpe')
box off; 

subplot(143)
bar(1,R.B,'FaceColor',color,'EdgeColor',color); hold on
line([0 2],[R.bounds.Ub.B R.bounds.Ub.B],'Color','k','LineWidth',1.5); hold on; %Ub
line([0 2],[R.bounds.Lb.B R.bounds.Lb.B],'Color','k','LineWidth',1.5); hold on; %Lb
title('B')
box off; 

subplot(144)
bar(1,R.kR,'FaceColor',color,'EdgeColor',color); hold on
%bar(2,R.kY,'FaceColor',color,'EdgeColor',color); hold on
line([0 2],[R.bounds.Ub.kR R.bounds.Ub.kR],'Color','k','LineWidth',1.5); hold on % Ub
line([0 2],[R.bounds.Lb.kR R.bounds.Lb.kR],'Color','k','LineWidth',1.5); hold on % Ub
title('kR')
box off 

sgtitle([info.subj ' Trial ', num2str(info.trial)])
end
