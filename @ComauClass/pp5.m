% [x, Dx, DDx] = pp5(t, xi, xf)
% [x, Dx, DDx] = pp5(t, xi, xf, Dxi, Dxf)
% [x, Dx, DDx] = pp5(t, xi, xf, Dxi, Dxf, DDxi, DDxf)
% ========================================================================
% x: dim-by-N array of position
% Dx: dim-by-N array of velocities
% DDx: dim-by-N array of accelerations
% t: time vector
% xi: dim-by-1 array of initial positions
% xf: dim-by-1 arrayo of final positions
% Dxi: dim-by-1 array of initial velocities     (optional)
% Dxf: dim-by-1 array of final velocities       (optional)
% DDxi: dim-by-1 array of initial accelerations (optional)
% DDxf: dim-by-1 array of final accelerations   (optional)
function [x, Dx, DDx] = pp5(t, xi, xf, varargin)
    t0 = t(1);
    t1 = t(end);
    A = [
        t0^5, t0^4, t0^3, t0^2, t0, 1;
        5*t0^4, 4*t0^3, 3*t0^2, 2*t0, 1, 0;
        20*t0^3, 12*t0^2, 6*t0, 2, 0, 0;
        t1^5, t1^4, t1^3, t1^2, t1, 1;
        5*t1^4, 4*t1^3, 3*t1^2, 2*t1, 1, 0;
        20*t1^3, 12*t1^2, 6*t1, 2, 0, 0];
    dim = numel(xi);
    N = numel(t);
    x = zeros(dim, N);
    Dx = zeros(dim, N);
    DDx = zeros(dim, N);
    if nargin == 3
        Dxi = zeros(dim, 1);
        Dxf = zeros(dim, 1);
        DDxi = zeros(dim, 1);
        DDxf = zeros(dim, 1);
    elseif nargin == 5
        Dxi = varargin{1};
        Dxf = varargin{2};
        DDxi = zeros(dim, 1);
        DDxf = zeros(dim, 1);
    elseif nargin == 7
        Dxi = varargin{1};
        Dxf = varargin{2};
        DDxi = varargin{3};
        DDxf = varargin{4};
    else
        error("Wrong arguments");
    end

    for i = 1:dim
        b = [xi(i); Dxi(i); DDxi(i); xf(i); Dxf(i); DDxf(i)];
        c = A\b;
        Dc = polyder(c);
        DDc = polyder(Dc);
        x(i,:) = polyval(c, t);
        Dx(i,:) = polyval(Dc, t);
        DDx(i,:) = polyval(DDc, t);
    end
end