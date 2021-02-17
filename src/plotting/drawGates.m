function [] = drawGates(gates,camera,measurement)
    nn = size(gates,1);
    % Group the projected points by their corresponding measurement:
    c = rand(nn,3);
    cla;
    f = camera.K(1,1);
    resx = camera.resolution(1);
    resy = camera.resolution(2);
    scatter(nan,nan,60,[0,0,0],'.'); hold on; axis equal; grid on
    scatter(nan,nan,60,[0,0,0],'x','LineWidth',2);
    plot(nan,nan,'-k');
    scatter(measurement(1:nn)*resx/f,measurement(1+nn:2*nn)*resy/f,60,c,'.');
    scatter(gates(:,1)*resx/f,gates(:,2)*resy/f,60,c,'x','LineWidth',2);
    for ii = 1:nn
        % Calculate the semi-major and semi-minor axes:
        x0 = gates(ii,1)*resx/f;
        y0 = gates(ii,2)*resy/f;
        sig3x = gates(ii,3)*resx/f;
        sig3y = gates(ii,4)*resy/f;
        top    = y0+sig3y;
        bottom = y0-sig3y;
        left   = x0-sig3x;
        right  = x0+sig3x;
        a = abs(top-bottom);
        b = abs(left-right);
        rotation = eye(2);
        drawEllipse(x0,y0,a,b,rotation,c(ii,:))
    end
    
    % Set the resolution of the image:
    legend('Measurement','Prediction','3-\sigma Gate')
    axis equal
    xlim([-resy/2 resy/2])
    ylim([-resy/2 resy/2])
%     axis equal
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
end