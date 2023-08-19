# rationalspline.m
A quick MATLAB snippet to implement shape-preserving rational spline Hermite interpolation.[(Cai and Judd, Economic Letters 2012)](https://www.sciencedirect.com/science/article/pii/S0165176512002558)
Requires MATLAB 2015a or later. (dependency: `discretize`)

## Example
```
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
```
<img width="618" alt="Screenshot 2023-08-19 at 7 24 01 AM" src="https://github.com/sunhamkim/rationalspline/assets/50336173/cef69bec-295f-41b6-a939-4d827ccf426f">
