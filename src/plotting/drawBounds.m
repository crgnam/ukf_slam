function [] = drawBounds(t,sig,num)
    t = t(1,1:num:end);
    sig = sig(1:num:end);
    FA = 0.2;
    a = area(t,[sig; 2*sig; 3*sig]'); hold on
    a(1).FaceColor = 'g';
    a(1).FaceAlpha = FA;
    a(1).EdgeColor = 'none';
    a(2).FaceColor = 'y';
    a(2).FaceAlpha = FA;
    a(2).EdgeColor = 'none';
    a(3).FaceColor = 'r';
    a(3).FaceAlpha = FA;
    a(3).EdgeColor = 'none';
    a = area(t,-[sig; 2*sig; 3*sig]');
    a(1).FaceColor = 'g';
    a(1).FaceAlpha = FA;
    a(1).EdgeColor = 'none';
    a(2).FaceColor = 'y';
    a(2).FaceAlpha = FA;
    a(2).EdgeColor = 'none';
    a(3).FaceColor = 'r';
    a(3).FaceAlpha = FA;
    a(3).EdgeColor = 'none';
end