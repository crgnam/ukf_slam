%% Show Map Initialization:
% Snippet of code to throw into orbitalA.m to get an animation of the map
% initialization projection step

figure('units','normalized','outerposition',[0 0 1 1])
orex.drawRays(image_lmks,radius_hat,r_hat);
bennu.drawBody();
% bennu.lght = [];
set(bennu.ptch,'AmbientStrength',.8,'FaceAlpha',0.5)
bennu.drawLmks(visible,'MarkerSize',30);
setLims(apoapsis)
camva(2)
legend('Rays','Approximate Sphere','Initialized Map','Bennu','Visible LMKs','Unseen LMKs','location','northeast')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
vid = VideoWriter('results/map_initialization.mp4','MPEG-4');
vid.Quality = 90;
open(vid)
for ii = 1:230
    view([30+ii/2 20])
    drawnow
    frame = getframe(gcf);
    writeVideo(vid,frame);
end
for ii = 1:230
    view([30+230/2-ii/2 20])
    drawnow
    frame = getframe(gcf);
    writeVideo(vid,frame);
end
close(vid)

