% Example and visualization

x = linspace(0.1,3,30).';   z = linspace(0.1,3,300).';
v = log(x);
s = 1./x;
[f,df,d2f] = rationalspline(x,z,v,s);

truf = log(z);
trudf = 1./z;

figure;
tt = tiledlayout(1,2);

nexttile
hold on;
plot(z,f);
plot(z,truf);
legend('Approximated Level','True Level','Location','best');
hold off;

nexttile
hold on;
plot(z,df);
plot(z,trudf);
legend('Approximated Slope','True Slope','Location','best');
hold off;

title(tt,'Interpolation of y = log(x)');
