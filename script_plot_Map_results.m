% (c) 2012 M Schaub -- michael.schaub09@imperial.ac.uk
variation =(L_exp-h_exp)./(h_exp);
%variation = [0 diff(L_exp)];
figure

subplot(2,1,1)
hold all
[ax, h1, h2]=plotyy(time,N_new,time,variation);
line(1,N,'Color','b','Marker','+','Parent',ax(1))
line(1,(L-h)/h,'Color',[0 127/255 0],'Marker','+','Parent',ax(2))
xlabel('Markov time');


set(ax(1),'YTickMode','auto','YTickLabelMode','auto','YAxisLocation','left');
set(get(ax(1),'Ylabel'),'String','# communities c');
set(get(ax(2),'Ylabel'),'String','compression gap \delta');

set(ax(1),'XLim', [10^floor(log10(time(1))) 10^ceil(log10(time(end)))], 'YLim', [0 max(N_new)*1.1], 'XScale','log');
set(ax(2),'XLim', [10^floor(log10(time(1))) 10^ceil(log10(time(end)))], 'YLim', [0 max([variation (L-h)/L])*1.1], 'XScale','log');
