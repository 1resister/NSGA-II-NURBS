function cntlpts = gencntlpts(p, Q, weightvec, U, temp)
% 兼容旧接口的控制点生成函数。
%
% 新项目中主流程使用 nurbs_global_smoothing.m。保留本函数是为了兼容旧
% NSGA-II/NURBS 调用入口：输入离散点 Q 后，按弦长参数自动生成控制点。

if size(Q, 1) ~= 3
    Q = Q';
end
points = Q';
nctrl = numel(weightvec);
if nargin >= 5 && ~isempty(temp)
    nctrl = temp;
end

s = [0; cumsum(vecnorm(diff(points), 2, 2))];
if s(end) < eps
    s = (0:size(points, 1)-1)';
end
target_s = linspace(0, s(end), nctrl);
cntlpts = interp1(s, points, target_s, 'pchip')';

if nargin >= 1 && p >= nctrl
    error('NURBS 阶数必须小于控制点数量。');
end
if nargin >= 4 && ~isempty(U)
    expected_knots = nctrl + p + 1;
    if numel(U) ~= expected_knots
        warning('节点向量长度与控制点/阶数不匹配，主流程将自动生成节点。');
    end
end
end
