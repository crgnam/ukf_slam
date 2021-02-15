% v = VideoReader('../results/media4.mp4');
% 
% vid = VideoWriter('../results/map_discovery2.mp4','MPEG-4');
% vid.Quality = 90;
% open(vid)
% 
% iter = 1;
% while hasFrame(v)
%     frame = readFrame(v);
%     frame2 = frame(50:end-100,300:end-150,:);
%     if mod(iter,10) == 0
%         writeVideo(vid,frame2);
%     end
%     iter = iter+1;
% end
% close(vid)
% close(v)

v = VideoReader('../results/map_discovery2.mp4');

vid = VideoWriter('../results/map_discovery3.mp4','MPEG-4');
vid.Quality = 50;
open(vid)

final = 8640;
iter = 1;
while hasFrame(v)
    frame = readFrame(v);
    writeVideo(vid,frame);
    iter = iter+1;
    disp(iter/final)
end
close(vid)
close(v)