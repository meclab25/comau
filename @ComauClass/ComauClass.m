classdef ComauClass < handle
    %% PUBLIC PROPERTIES
    properties (Access = public)
        q0   (3, 1) double;     % Joint variables
        r0   (3, 1) double;     % Position of end-effector
        phi0 (3, 1) double;     % Euler angles of end-effector
        % T = Tf(q)
        % =================================================================
        Tf  function_handle;    % Function handle
        % -----------------------------------------------------------------
        % V = Vf(q, Dq)
        % =================================================================
        Vf  function_handle;    % Function hangle for twist
        % -----------------------------------------------------------------
        % DV = DVf(q, Dq, DDq)
        % =================================================================
        DVf function_handle;    % Function handle for rate of twist
        % -----------------------------------------------------------------
        % Dq = Dqf(q, r, Dr)
        % =================================================================
        Dqf function_handle;    % Function handle for inverse velocity
        % -----------------------------------------------------------------
        % DDq = DDqf(q, Dq, r, Dr, DDr)
        % =================================================================
        DDqf function_handle;   % Function handle for inverse acceleration
        % -----------------------------------------------------------------

    end
    %% PUBLIC PROPERTIES - SIMULATION ARRAYS
    properties (Access = public)
        t   double; % time array
        q   double; % joint variables
        Dq  double; % joint velocities
        DDq double; % joint accelerations
        r   double; % end-effector position
        Dr  double; % end-effector velocity
        DDr double; % end-effector acceleration
        phi double; % end-effector euler angles
        w   double; % end-effector angular velocity
        Dw  double; % end-effector angular acceleration
    end
    %% CONSTANT PROPERTIES
    properties (Constant = true)
        L1 = 350e-3;
        L2 = 830e-3;
        L3 = 1160e-3;
        L4 = 2269.18e-3;
        % L4 = 1492.18e-3;
        L5 = 250e-3;
        pX = 2.40;
        pY = 0.35;
    end
    %% HIDDEN PROPERTIES
    properties (Hidden = true)
        Tj2f function_handle;
        Tj3f function_handle;
        TPf function_handle;
        fig
        an1
        an2
        an3
        an4
        an5
        an6
        an7
        an8

    end
    %% HIDDEN/CONSTANT PROPERTIES
    properties (Constant = true, Hidden = true)
        x = [1,0,0].';
        y = [0,1,0].';
        z = [0,0,1].';
    end
    %% PUBLIC METHODS
    methods (Access = public)
        % obj = ComauClass(q0)
        % =================================================================
        % obj: class object
        % q0: initial joint angles
        % -----------------------------------------------------------------
        % Constructor of the class
        function obj = ComauClass(q0)
            q = sym('q', [3, 1], 'real');
            Dq = sym('Dq', [3, 1], 'real');
            DDq = sym('DDq', [3, 1], 'real');
            % Forward Kinematics
            r = [obj.L1 - obj.L5; 0; obj.L2 + obj.L3 + obj.L4];
            T0 = [obj.Rot(obj.y, -pi/2), r; 0, 0, 0, 1];
            S1 = [-obj.z; zeros(3,1)];
            S2 = [obj.y; obj.so3(-obj.y)*(obj.L1*obj.x + obj.L2*obj.z)];
            S3 = [-obj.y; obj.so3(obj.y)*(obj.L1*obj.x + (obj.L2 + obj.L3)*obj.z)];
            T = obj.mSE3(S1, q(1)) * obj.mSE3(S2, q(2)) * obj.mSE3(S3, q(2)+q(3)) * T0;
            obj.Tf = matlabFunction(T, 'vars', {q});
            % Veloctiy
            T1 = obj.mSE3(S1, q(1));
            T2 = T1*obj.mSE3(S2, q(2));
            % T3 = T2*obj.mSE3(S3, q(2)+q(3));
            Jac = [S1, obj.ad3(T1)*S2, obj.ad3(T2)*S3];
            V = Jac*(Dq + [0; 0; Dq(2)]);
            obj.Vf = matlabFunction(V, 'vars', {q, Dq});
            obj.q0 = q0;
            T0 = obj.Tf(q0);
            obj.r0 = T0(1:3, 4);
            obj.phi0 = obj.EulXYZ(T0(1:3, 1:3));
            % Acceleration
            DV = jacobian(V, [q; Dq])*[Dq; DDq];
            obj.DVf = matlabFunction(DV, 'vars', {q, Dq, DDq});

            % Inverse velocity and acceleration
            r = sym('r', [3, 1], 'real');
            Dr = sym('Dr', [3, 1], 'real');
            DDr = sym('DDr', [3, 1], 'real');

            eq1 = V(4:end) + obj.so3(V(1:3))*r - Dr;
            Dqs = solve(eq1, Dq);
            eq2 = DV(4:6) + obj.so3(DV(1:3))*r + obj.so3(V(1:3))*Dr - DDr;
            DDqs = solve(eq2, DDq);

            Dqf = [Dqs.Dq1; Dqs.Dq2; Dqs.Dq3];
            DDqf = [DDqs.DDq1; DDqs.DDq2; DDqs.DDq3];
            obj.Dqf = matlabFunction(Dqf, 'vars', {q, r, Dr});
            obj.DDqf = matlabFunction(DDqf, 'vars', {q, Dq, r, Dr, DDr});

            % Visualization for position of each joint
            Tj20 = [eye(3), [obj.L1, 0, obj.L2].'; 0, 0, 0, 1];
            Tj2 = obj.mSE3(S1, q(1))*Tj20;
            Tj30 = Tj20 * [eye(3), [0, 0, obj.L3].'; 0, 0, 0, 1];
            Tj3 = obj.mSE3(S1, q(1))*obj.mSE3(S2, q(2))*Tj30;
            TP0 = [eye(3), [obj.L1 - obj.L5, 0, obj.L2 + obj.L3].'; 0, 0, 0, 1];
            TP = obj.mSE3(S1, q(1))*obj.mSE3(S2, q(2))*obj.mSE3(S3, q(2)+q(3))*TP0;
            obj.Tj2f = matlabFunction(Tj2, 'vars', {q});
            obj.Tj3f = matlabFunction(Tj3, 'vars', {q});
            obj.TPf = matlabFunction(TP, 'vars', {q});
        end
        % obj = Animate
        % =================================================================
        % -----------------------------------------------------------------
        function Animate(obj, varargin)
            if nargin == 2
                r0 = varargin{1};
            else
                r0 = [0,0,0].';
            end
            Tj20 = obj.Tj2f(obj.q0);
            Tj30 = obj.Tj3f(obj.q0);
            TP0 = obj.TPf(obj.q0);
            T0 = obj.Tf(obj.q0);
            p01 = r0 + Tj20(1:3, 4);
            p02 = r0 + Tj30(1:3, 4);
            p03 = r0 + TP0(1:3, 4);
            p04 = r0 + T0(1:3, 4);
            obj.fig = figure(1);
            obj.an1 = animatedline([r0(1), p01(1)], ...
                [r0(2), p01(2)], ...
                [r0(3), p01(3)], ...
                "Color", "k", "LineWidth", 2);
            obj.an2 = animatedline([p01(1), p02(1)], ...
                [p01(2), p02(2)], ...
                [p01(3), p02(3)], ...
                "Color", "m", "LineWidth", 2);
            obj.an3 = animatedline([p02(1), p03(1)], ...
                [p02(2), p03(2)], ...
                [p02(3), p03(3)], ...
                "Color", "k", "LineWidth", 2);
            obj.an4 = animatedline([p03(1), p04(1)], ...
                [p03(2), p04(2)], ...
                [p03(3), p04(3)], ...
                "Color", "m", "LineWidth", 2);
            obj.an5 = animatedline(p04(1), p04(2), p04(3), ...
                "Color", "r", "LineWidth", 2);
            c = 0.5;
            obj.an6 = animatedline([p04(1), p04(1)+c*T0(1, 1)], ...
                [p04(2), p04(2)+c*T0(2, 1)], ...
                [p04(3), p04(3)+c*T0(3, 1)], ...
                "Color", "r", "LineWidth", 1.2);
            obj.an7 = animatedline([p04(1), p04(1)+c*T0(1, 2)], ...
                [p04(2), p04(2)+c*T0(2, 2)], ...
                [p04(3), p04(3)+c*T0(3, 2)], ...
                "Color", "g", "LineWidth", 1.2);
            obj.an8 = animatedline([p04(1), p04(1)+c*T0(1, 3)], ...
                [p04(2), p04(2)+c*T0(2, 3)], ...
                [p04(3), p04(3)+c*T0(3, 3)], ...
                "Color", "b", "LineWidth", 1.2);

            % add patch
            temp.x = 2.40;
            temp.y = 0.35;
            temp.v1 = [temp.x; -temp.y; 0];
            temp.v2 = [temp.x; temp.y; 0];
            temp.v3 = obj.Rot(obj.z, 120*pi/180)*temp.v1;
            temp.v4 = obj.Rot(obj.z, 120*pi/180)*temp.v2;
            temp.v5 = obj.Rot(obj.z, 240*pi/180)*temp.v1;
            temp.v6 = obj.Rot(obj.z, 240*pi/180)*temp.v2;
            temp.P = [temp.v1, temp.v2, temp.v3, temp.v4, temp.v5, temp.v6];

            temp.T = [eye(3), [-1.2;0;0]; 0, 0, 0, 1];
            temp.P = temp.T*[temp.P; ones(1, 6)];
            temp.P = temp.P(1:3, :);
            patch(temp.P(1,:), temp.P(2,:), temp.P(3,:), 'cyan', 'FaceAlpha', 0.5)
            line([0, 0.5], [0, 0], [0, 0], 'Color', 'r', 'LineWidth', 2);
            line([0, 0], [0, 0.5], [0, 0], 'Color', 'g', 'LineWidth', 2);
            line([0, 0], [0, 0], [0, 0.5], 'Color', 'b', 'LineWidth', 2);
            view(120, 20); grid on; axis equal;
            axis([-3, 3, -3, 3, -2, 4]);
        end
        % PlotMe(ind, clearFlag)
        % =================================================================
        % ind: integer - index of trajectory
        % clearFlag: boolean - if true, then clear the line on update
        % -----------------------------------------------------------------
        function plotUpdate(obj, N, clearFlag, varargin)
            if N > numel(obj.t)
                warning("N to large");
            else
                if nargin == 4
                    r = varargin{1};
                else
                    r = [0,0,0].';
                end
                q_ = obj.q(:, N);
                Tj2 = obj.Tj2f(q_);
                Tj3 = obj.Tj3f(q_);
                TP = obj.TPf(q_);
                T = obj.Tf(q_);
                p01 = r + Tj2(1:3, 4);
                p02 = r + Tj3(1:3, 4);
                p03 = r + TP(1:3, 4);
                p04 = r + T(1:3, 4);
                if clearFlag
                    clearpoints(obj.an1);
                    clearpoints(obj.an2);
                    clearpoints(obj.an3);
                    clearpoints(obj.an4);
                    clearpoints(obj.an6);
                    clearpoints(obj.an7);
                    clearpoints(obj.an8);
                end
                c = 0.5;
                addpoints(obj.an1, [r(1), p01(1)], [r(2), p01(2)], [r(3), p01(3)])
                addpoints(obj.an2, [p01(1), p02(1)], [p01(2), p02(2)], [p01(3), p02(3)])
                addpoints(obj.an3, [p02(1), p03(1)], [p02(2), p03(2)], [p02(3), p03(3)])
                addpoints(obj.an4, [p03(1), p04(1)], [p03(2), p04(2)], [p03(3), p04(3)])
                addpoints(obj.an5, p04(1), p04(2), p04(3));
                addpoints(obj.an6, [p04(1), p04(1)+c*T(1, 1)], [p04(2), p04(2)+c*T(2, 1)], ...
                    [p04(3), p04(3)+c*T(3, 1)]);
                addpoints(obj.an7, [p04(1), p04(1)+c*T(1, 2)], [p04(2), p04(2)+c*T(2, 2)], ...
                    [p04(3), p04(3)+c*T(3, 2)]);    
                addpoints(obj.an8, [p04(1), p04(1)+c*T(1, 3)], [p04(2), p04(2)+c*T(2, 3)], ...
                    [p04(3), p04(3)+c*T(3, 3)]);
                drawnow;
            end
        end
        function plotMe(obj)
            obj.fig;
            plot3(obj.r(1,:), obj.r(2,:), obj.r(3,:));
        end
        % obj = FKTraj(time, qf, type)
        % =================================================================
        % obj: class object
        % time: time vector
        % qf: final joint angle
        % type: 'pp3', 'pp5', 'pp7', 'bcb'
        % -----------------------------------------------------------------
        % Forward Kinematics to populate a trajectory
        function obj = FKTraj(obj, time, qf, type)
            N = numel(time);
            r_ = zeros(3, N);
            Dr_ = zeros(3, N);
            DDr_ = zeros(3, N);
            phi_ = zeros(3, N);
            w_ = zeros(3, N);
            Dw_ = zeros(3, N);
            if strcmp(type, 'pp3')
                [q_, Dq_, DDq_] = obj.pp3(time, obj.q0, qf);
            elseif strcmp(type, 'pp5')
                [q_, Dq_, DDq_] = obj.pp5(time, obj.q0, qf);
            elseif strcmp(type, 'pp7')
                [q_, Dq_, DDq_] = obj.pp7(time, obj.q0, qf);
            elseif strcmp(type, 'bcb')
                warning("Splitting time into arbitrary pieces")
                t1 = time(end)/4;
                t2 = t1;
                t3 = t1;
                [q_, Dq_, DDq_] = obj.bcb(time, t1, t2, t3, obj.q0, qf);
            else
                error("Wrong type for FKTraj")
            end
            for i = 1:N
                T = obj.Tf(q_(:, i));
                V = obj.Vf(q_(:, i), Dq_(:, i));
                DV = obj.DVf(q_(:, i), Dq_(:, i), DDq_(:, i));
                r_(:, i) = T(1:3, 4);
                Dr_(:, i) = V(4:6) + obj.so3(V(1:3))*r_(:, i);
                DDr_(:, i) = DV(4:6) + obj.so3(DV(1:3))*r_(:, i) + ...
                    obj.so3(V(1:3))*Dr_(:, i);
                phi_(:, i) = obj.EulXYZ(T(1:3, 1:3));
                w_(:, i) = V(1:3);
                Dw_(:, i) = DV(1:3);
            end
            obj.t = time;
            obj.q = q_;
            obj.Dq = Dq_;
            obj.DDq = DDq_;
            obj.r = r_;
            obj.Dr = Dr_;
            obj.DDr = DDr_;
            obj.phi = phi_;
            obj.w = w_;
            obj.Dw = Dw_;
        end

        % obj = IKTraj(time, p, Dp, DDp)
        % =================================================================
        % time: time vector
        % p: 3-by-N array of displacement (delta r)
        % Dp: 3-by-N array of velocity
        % DDp: 3-by-N array of acceleration
        function obj = IKTraj(obj, time, p, Dp, DDp)
            N = numel(time);
            r_ = p;
            Dr_ = Dp;
            DDr_ = DDp;
            q_ = zeros(3, N);
            Dq_ = zeros(3, N);
            DDq_ = zeros(3, N);
            phi_ = zeros(3, N);
            w_ = zeros(3, N);
            Dw_ = zeros(3, N);
            for i = 1:N
                r__ = r_(:, i);
                Dr__ = Dr_(:, i);
                DDr__ = DDr_(:, i);
                q__ = obj.invKin(r__);
                Dq__ = obj.Dqf(q__, r__, Dr__);
                DDq__ = obj.DDqf(q__, Dq__, r__, Dr__, DDr__);
                V = obj.Vf(q__, Dq__);
                DV = obj.DVf(q__, Dq__, DDq__);
                T = obj.Tf(q__);
                phi_(:, i) = obj.EulXYZ(T(1:3, 1:3));
                q_(:, i) = q__;
                Dq_(:, i) = Dq__;
                DDq_(:, i) = DDq__;
                w_(:, i) = V(1:3);
                Dw_(:, i) = DV(1:3);
            end
            obj.t = time;
            obj.q = q_;
            obj.Dq = Dq_;
            obj.DDq = DDq_;
            obj.r = r_;
            obj.Dr = Dr_;
            obj.DDr = DDr_;
            obj.phi = phi_;
            obj.w = w_;
            obj.Dw = Dw_;
        end
    end

    methods (Access = private)
        % phi = EulXYZ(R)
        % =================================================================
        % phi: 3-by-1 vector of Euler angles - XYZ order
        % R: 3-by-3 matrix in SO(3) -  XYZ order
        % -----------------------------------------------------------------
        function phi = EulXYZ(obj, R)
            pz = atan2(-R(1,2), R(1,1));
            px = atan2(-R(2,3), R(3,3));
            cpsi = cos(pz);
            cth = R(1,1)/cpsi;
            sth = R(1,3);
            py = atan2(sth, cth);
            phi = [px; py; pz];
        end
        % J = ad3(T)
        % =================================================================
        % J: 6-by-6 matrix
        % T: 4-by-4 matrix in SE(3)
        % -----------------------------------------------------------------
        % Adjoint map between a twist V and joint velocity Dq
        function J = ad3(obj, T)
            R = T(1:3, 1:3);
            p = T(1:3, 4);
            pS = obj.so3(p);
            J = [R, zeros(3,3); pS*R, R];
        end
        % R = Rot(w, theta)
        % =================================================================
        % R: 3-by-3 matrix in SO(3)
        % w: 3-by-1 unit vector
        % theta: scalar (angle of rotation)
        % -----------------------------------------------------------------
        % Rotation matrix for unit vector w and angle theta
        function R = Rot(obj, w, theta)
            if norm(w) ~= 1
                warning("non-unit vector");
                w = w/norm(w);
            end
            wS = obj.so3(w);
            I = eye(3);
            s = sin(theta);
            c = cos(theta);
            R = I + s*wS + (1 - c)*wS^2;
        end
        % vS = so3(v)
        % =================================================================
        % vS: 3-by-3 matrix in so(3)
        % v: 3-by-1 vector
        % -----------------------------------------------------------------
        % Skew symmetric matrix form of a vector.
        % Enables cross product (v x u) = vS*u
        function vS = so3(~, v)
            vS = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
        end
        % SS = se3(S)
        % =================================================================
        % S: 6-by-1 screw
        % SS: 4-by-4 matrix in se(3)
        % -----------------------------------------------------------------
        % S = [w; v], where |w| = 1 and v = -w x p.
        % w is an angle of rotation and p is the vector from {s} to w.
        % Map a screw to the argument of the matrix exponential in SE(3).
        function SS = se3(obj, S)
            w_ = S(1:3);
            v = S(4:6);
            SS = [obj.so3(w_), v; 0, 0, 0, 0];
        end
        % T = mSE3(S, theta)
        % =================================================================
        % T: 4-by-4 matrix in SE(3)
        % S: 6-by-1 screw vector
        % theta: angle of rotation about the screw
        % -----------------------------------------------------------------
        % Maps a screw to SE(3) using the matrix exponential.
        function T = mSE3(obj, S, theta)
            w_ = S(1:3);
            v = S(4:6);
            if norm(w_) ~= 1
                warning("non-unit vector")
                w_ = w_/norm(w_);
            end
            R = obj.Rot(w_, theta);
            I = eye(3);
            s = sin(theta);
            c = cos(theta);
            wS = obj.so3(w_);
            r_ = (I*theta + (1 - c)*wS + (theta - s)*wS^2)*v;
            T = [R, r_; 0,0,0,1];
        end
        % q = invKin(r)
        % =================================================================
        % q: 3-by-1 vector of joint angles
        % r: 3-by-1 vector of end-effector position values
        % -----------------------------------------------------------------
        function q = invKin(obj, r)
            q1 = -atan2(r(2), r(1));
            r1 = obj.Rot(obj.z, q1)*r;
            r2 = r1 - [obj.L1; 0; obj.L2];
            x_ = r2(1);
            z_ = r2(3);
            L6 = sqrt(obj.L4^2 + obj.L5^2);
            psi = atan(obj.L5/obj.L4);
            L = sqrt(x_^2 + z_^2);
            theta = atan(z_/x_);
            phi6 = acos((obj.L3^2 + L^2 - L6^2)/(2*obj.L3*L));
            q2 = pi/2 - theta - phi6;
            phi_ = acos((L6^2+ obj.L3^2 - L^2) / (2*L6*obj.L3));
            q3star = pi + psi - phi_;
            q3 = -q3star - q2;
            q = [q1; q2; q3];
            if any(imag(q) ~= 0)
                fprintf("Singularity:\n\t x = %.2f\t y = %.2f\t z = %.2f\n", r);
                error("Singularity");
            end
        end

    end
    %% STATIC METHODS
    methods (Static)
        % [x, Dx, DDx] = pp3(t, xi, xf)
        % [x, Dx, DDx] = pp3(t, xi, xf, Dxi, Dxf)
        % =================================================================
        % x: dim-by-N array of position
        % Dx: dim-by-N array of velocities
        % DDx: dim-by-N array of accelerations
        % t: time vector
        % xi: dim-by-1 array of initial positions
        % xf: dim-by-1 arrayo of final positions
        % Dxi: dim-by-1 array of initial velocities (optional)
        % Dxf: dim-by-1 array of final velocities   (optional)
        % -----------------------------------------------------------------
        % 3-rd order polynomial trajectory
        [x, Dx, DDx] = pp3(t, xi, xf, varargin);

        % [x, Dx, DDx] = pp5(t, xi, xf)
        % [x, Dx, DDx] = pp5(t, xi, xf, Dxi, Dxf)
        % [x, Dx, DDx] = pp5(t, xi, xf, Dxi, Dxf, DDxi, DDxf)
        % =================================================================
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
        % -----------------------------------------------------------------
        % 5'th order polynomial trajectory generator
        [x, Dx, DDx] = pp5(t, xi, xf, varargin);

        % [x, Dx, DDx] = pp7(t, xi, xf)
        % [x, Dx, DDx] = pp7(t, xi, xf, Dxi, Dxf)
        % [x, Dx, DDx] = pp7(t, xi, xf, Dxi, Dxf, DDxi, DDxf)
        % [x, Dx, DDx] = pp7(t, xi, xf, Dxi, Dxf, DDxi, DDxf, DDDxi, DDDxf)
        % =================================================================
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
        % DDDxi: dim-by-1 array of initial jerk         (optional)
        % DDDxf: dim-by-1 array of final jerk           (optional)
        % -----------------------------------------------------------------
        % 7'th order polynomial trajectory generator
        [x, Dx, DDx] = pp7(t, xi, xf, varargin);

        % [x, Dx, DDx] = bcb(t, t1, t2, t3, xi, xf)
        % =================================================================
        % x: dim-by-N array of position
        % Dx: dim-by-N array of velocities
        % DDx: dim-by-N array of accelerations
        % t: time vector
        % t1: bang-time (ramping velocity from zero)
        % t2: coast-time (constant velocity)
        % t3: bang-time (ramping velocity to zero)
        % xi: dim-by-1 array of initial position
        % xf: dim-by-1 array of final position
        % -----------------------------------------------------------------
        % Bang-coast-bang: start and end at zero velocity
        [x, Dx, DDx] = bcb(t, t1, t2, t3, xi, xf);

    end
end