% visualization parameters

set(gcf,'color','w')
set(gcf,'paperunits','inches') 
set(gcf,'units','inches') 
set(gcf,'PaperPositionMode', 'manual');
set(gcf,'papersize',[widthIn+marginleft+marginright,heightIn+marginbottom+margintop])
set(gcf,'paperposition',[0,0,widthIn+marginleft+marginright,heightIn+marginbottom+margintop])
set(gcf,'renderer', 'painters');
set(gcf,'Position',[marginleft,marginbottom,widthIn,heightIn],'units','inches') 