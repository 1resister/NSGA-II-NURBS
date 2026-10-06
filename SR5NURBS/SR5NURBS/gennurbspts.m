function nurbspts = gennurbspts(cntlpts, curdeg, weightsvec, knotvec, N)
% 兼容旧接口的 NURBS 曲线采样函数。
%
% 新主流程使用 evaluate_nurbs_curve.m，本函数仅负责把旧的 3 x n 控制点
% 格式转换为新实现需要的 n x 3 格式。

if size(cntlpts, 1) == 3
    control_points = cntlpts';
else
    control_points = cntlpts;
end

u = linspace(0, 1, N);
points = evaluate_nurbs_curve(control_points, weightsvec(:)', curdeg, knotvec(:)', u);
nurbspts = points';
end
