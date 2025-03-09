clear; close all; clc;
Ts = 1e-3;
time = 0:Ts:5;
q0 = [0, 0, -pi].';
cm = ComauClass(q0);
cm = cm.FKTraj(time, [pi, 0, -pi/2.].', 'pp7');
cm.Animate;
for i = 1:5:5000
    cm.plotUpdate(i, 1)
end

%% Inverse Kinematics
rf = [0; -2; 3];
ri = cm.r(:, end);
q0 = cm.q(:, end);
[r, Dr, DDr] = cm.pp5(time, ri, rf);
cm = ComauClass(q0);
cm = cm.IKTraj(time, r, Dr, DDr);
cm.Animate;
for i = 1:50:5000
    cm.plotUpdate(i, 1)
end



