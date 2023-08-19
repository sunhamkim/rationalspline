% x, v, s: m x 1 column vector
% output.coef: mx4 matrix
% output.interval: x = mx1 vector

function [f, df, d2f] = rationalspline(x,z,v,s)

% MATLAB Implmentation of shape-preserving rational spline Hermite
% interpolation, proposed by:
% Cai, Y., & Judd, K. L. (2012). 
% Dynamic programming with shape-preserving rational spline Hermite interpolation. 
% Economics Letters, 117(1), 161-164.
%
% Author: Sunham Kim
%
% Input:  A function v and its first derivative s, on a sample grid x.
% Output: An interpolant f and its first and second derivatives df and d2f
%         on a full grid z.
%
% Usage:
%   [f, df, d2f] = rationalspline(x,z,v,s);
% 
% 
% Example:
%   x = linspace(0.1,3,20);   z = linspace(0.1,3,200)
%   v = log(x);
%   s = 1./x;
%   [f,df,d2f] = defrspline(x,z,v,s);
%
% Usage restrictions
%    x, v, s are kx1 column vectors and supposed to be sorted with respect
%    to x.
%    v must be concave and monotonically increasing
%    extrapolation is not supported; min(x) <= min(z) and max(z) <= max(x)
%



% Dimension and Boundary Check
if any([~iscolumn(x),~iscolumn(z),~iscolumn(v),~iscolumn(s)]), error('All inputs must be a column vector.'); end
if min(x(:)) > min(z(:)), error('min(z) must be greater than or equal to min(x).'); end
if max(x(:)) < max(z(:)), error('max(z) must be smaller than or equal to max(x).'); end
if ~issorted(x,'ascend'), error('x has to be sorted in ascending order; all other arrays have to be sorted with respect to x.'); end


% Obtain Rational Spline Coefficients
c1 = v(1:end-1);
c2 = diff(v)./diff(x);
c3 = s(1:end-1) - c2;
c4 = s(2:end)   - c2;
coef = [c1,c2,c3,c4];

% Evaluation Step
ix = discretize(z(:),x(:));
zL = z-x(ix);               % x-x_i. Left.
zR = z-x(1+ix);             % x-x_i+1. Right.
coefGrid = coef(ix,:);      % c_i

f = coefGrid(:,1)     + ...
    coefGrid(:,2).*zL + ...
    (prod(coefGrid(:,3:4),2).*zL.*zR ) ./ ( coefGrid(:,3).*zL + coefGrid(:,4).*zR );

df = coefGrid(:,2) + ...
        prod(coefGrid(:,3:4),2).*( coefGrid(:,3).*(zL.^2) + coefGrid(:,4).*(zR.^2) ) ./ ...
            ( coefGrid(:,3).*zL + coefGrid(:,4).*zR ).^2;

d2f = 2*(prod(coefGrid(:,3:4),2).^2).*((zR-zL).^2) ./ ...
            (-( coefGrid(:,3).*zL + coefGrid(:,4).*zR ).^3);

end
