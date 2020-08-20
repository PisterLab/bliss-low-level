function [] = animateSolarSystem(xtrajs,ytrajs,ztrajs,tSim,recordAnimation2File)

lnwidth = 1.0;
fsize = 12;
mksize = 50;
% rcsize = mksize+50;
% rcalpha = 0.25;

sec_in_solar_day = 86400;

m_AU   = (149597870.7)*10^3;  % [m/AU] Length of an Astronomical Unit in meters
xtrajs = xtrajs/m_AU;         % normalize by AU
ytrajs = ytrajs/m_AU;         % normalize by AU
ztrajs = ztrajs/m_AU;         % normalize by AU

m    = size(xtrajs,1);
tmax = size(xtrajs,2);

if recordAnimation2File
    vfname = ['animation_' datestr(now,'mm_dd_yyyy_HH_MM_SS_AM')]; %#ok<*UNRCH>
    v = VideoWriter(fullfile(fileparts(which(mfilename)),'animation',[vfname '.mp4']),'MPEG-4');
    v.Quality = 80;
    v.FrameRate = 20;
    open(v);
end

grid on;
hold on;
axis equal;
xlabel('1-Dir [AU]'); ylabel('2-Dir [AU]'); zlabel('3-Dir [AU]');
set(gca, 'FontSize', fsize,'FontWeight','bold')
set(gcf,'color','w');
% view(-25.6472,15.4648);
view(0,90);
% view(-13.8411,13.5651);
xlimmin = min(xtrajs,[],'all');
xlimmax = max(xtrajs,[],'all');
ylimmin = min(ytrajs,[],'all');
ylimmax = max(ytrajs,[],'all');
xlim([xlimmin xlimmax]);
ylim([ylimmin ylimmax]);

scatter3(0,0,0,10*mksize,[243 130 53]/255,'filled'); % Sun

I  = cell(m,1);
Ih = cell(m,1);

for i = 1:m % using a loop just in case we decide to add more bodies later
    if i == 34
                I{i} = animatedline('linewidth',0.2*lnwidth,'color',[0 0.447 0.741]);
        Ih{i} = scatter3(nan,nan,nan,2*mksize,[0 0.447 0.741],'filled');    % Pale Blue Dot
    elseif i > 31 && i~=34
        I{i} = animatedline('linewidth',0.2*lnwidth,'color',[170 170 170]/255);
        Ih{i} = scatter3(nan,nan,nan,2*mksize,[170 170 170]/255,'filled');    % gray planets
    else
        I{i} = animatedline('linewidth',0.2*lnwidth);
        Ih{i} = scatter3(nan,nan,nan,0.5*mksize,'filled'); % NEOs
    end
end

step = 500; % No need to plot every single point. Takes too long.
for t = [1:step:tmax,tmax]
    days = text(-1.5,-3,0, [num2str((tSim(t)/sec_in_solar_day),'%0.1f') ' Solar Days']);
    days.FontSize = 14;
    days.HorizontalAlignment = 'center'; days.VerticalAlignment = 'top';
    
    for i = 1:m
        x = xtrajs(i,:);
        y = ytrajs(i,:);
        z = ztrajs(i,:);
        addpoints(I{i},x(t),y(t),z(t));
        set(Ih{i},'xdata',x(t),'ydata',y(t),'zdata',z(t))
%         if i == 2
%             if abs(t - idx_CPA) < step
%                 scatter3(x(t),y(t),z(t),rcsize,[0.850 0.325 0.0980],'filled','markerfacealpha',rcalpha) % Mark the CPA
%             end
%         end
    end
    
    drawnow limitrate;  % Too slow. Make the animation even faster.
    
    %Write Video
    if recordAnimation2File
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
    
    if t<tmax
        delete(days)
    end
    
end

if recordAnimation2File
    close(v);
end

end