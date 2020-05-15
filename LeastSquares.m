H = [160;155;143;162;170;175;167;163;156;172;180;159;169;157;161;171;181;
    177;140]
W = [68;78;56;70;72;78;68;54;43;71;65;59;63;62;90;45;62;57;80];

%w = ah^2 + bh + c

A = [H.^2, H, ones(length(H), 1)];
disp(A);
disp(transpose(A)*A)
disp(transpose(A)*W)
sol = inv(transpose(A)*A)*transpose(A)*W;

disp(sol);
a = sol(1, 1);
b = sol(2, 1);
c = sol(3, 1);
disp(a*185*185 + b*185 + c);

sol2 = pinv(A)*W;
disp(sol2);
