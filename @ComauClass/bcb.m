% [x, Dx, DDx] = bcb(t, t1, t2, t3, xi, xf)
% =======================================================================
% x: dim-by-N array of position
% Dx: dim-by-N array of velocities
% DDx: dim-by-N array of accelerations
% t: time vector
% t1: bang-time (ramping velocity from zero)
% t2: coast-time (constant velocity)
% t3: bang-time (ramping velocity to zero)
% xi: dim-by-1 array of initial position 
% xf: dim-by-1 array of final position
% -----------------------------------------------------------------------
% Bang-coast-bang: start and end at zero velocity
function [x, Dx, DDx] = bcb(t, t1, t2, t3, xi, xf)
    if t(end) >= t1 + t2 + t3
        M = [t1*(t1/2 + t2 + t3), t3^2/2; t1, t3];
        el1 = t < t1;
        el2 = t >= t1 & t < (t1 + t2);
        el3 = t >= (t1 +t2) & t < (t1 + t2 + t3);
        el4 = t >= (t1 + t2 + t3);
        dim = numel(xi);
        N = numel(t);
        x = zeros(dim, N);
        Dx = zeros(dim, N);
        DDx = zeros(dim, N);
        for i = 1:dim
            b = [xf(i) - xi(i); 0];
            c = M\b;
            c1 = c(1);
            c2 = c(2);

            x1 = c1*t(el1).^2/2;
            x2 = c1*t1^2/2 + c1*t1*(t(el2)-t1);
            x3 = c1*t1^2/2 + c1*t1*(t(el3)-t1) + c2*(t(el3)-(t1+t2)).^2/2;
            x4 = c1*t1*(t1/2+t2+t3) + c2*t3^2/2 + 0*t(el4);

            v1 = c1*t(el1);
            v2 = c1*t1 + 0*t(el2);
            v3 = c1*t1 + c2*(t(el3)- (t1+ t2));
            v4 = 0*t(el4);

            a1 = c1 + 0*t(el1);
            a2 = 0*t(el2);
            a3 = c2 + 0*t(el3);
            a4 = 0*t(el4);

            x(i,:) = xi(i) + [x1, x2, x3, x4];
            Dx(i,:) = [v1, v2, v3, v4];
            DDx(i,:) = [a1, a2, a3, a4];
        end

    else
        error("time-vector to short")
    end