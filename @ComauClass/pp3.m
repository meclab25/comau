% [x, Dx, DDx] = pp3(t, xi, xf)
% [x, Dx, DDx] = pp3(t, xi, xf, Dxi, Dxf)
% ========================================================================
% x: dim-by-N array of position
% Dx: dim-by-N array of velocities
% DDx: dim-by-N array of accelerations
% t: time vector
% xi: dim-by-1 array of initial positions
% xf: dim-by-1 arrayo of final positions
% Dxi: dim-by-1 array of initial velocities (optional)
% Dxf: dim-by-1 array of final velocities   (optional)
function [x, Dx, DDx] = pp3(t, xi, xf, varargin)
    t0 = t(1);
    t1 = t(end);
    A = [
        t0^3, t0^2, t0, 1;
        3*t0^2, 2*t0, 1, 0;
        t1^3, t1^2, t1, 1;
        3*t1^2, 2*t1, 1, 0];
    dim = numel(xi);
    N = numel(t);
    x = zeros(dim, N);
    Dx = zeros(dim, N);
    DDx = zeros(dim, N);
    if nargin == 3
        Dxi = zeros(dim, 1);
        Dxf = zeros(dim, 1);
    elseif nargin == 5
        Dxi = varargin{1};
        Dxf = varargin{2};
    else
        error("Wrong arguments")
    end

    for i = 1:dim
        b = [xi(i); Dxi(i); xf(i); Dxf(i)];
        c = A\b;
        Dc = polyder(c);
        DDc = polyder(Dc);
        x(i,:) = polyval(c, t);
        Dx(i,:) = polyval(Dc, t);
        DDx(i,:) = polyval(DDc, t);
    end
end