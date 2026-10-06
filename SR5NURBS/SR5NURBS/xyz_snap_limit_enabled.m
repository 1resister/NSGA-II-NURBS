function enabled = xyz_snap_limit_enabled(cfg)
%XYZ_SNAP_LIMIT_ENABLED Return whether the XYZ snap hard limit is enforced.
%
% Default to true for configurations created before the switch existed.

enabled = true;
if isfield(cfg, 'constraints') && ...
        isfield(cfg.constraints, 'enable_xyz_snap_limit')
    value = cfg.constraints.enable_xyz_snap_limit;
    if ~isscalar(value) || ...
            ~(islogical(value) || (isnumeric(value) && isfinite(value) && ...
            (value == 0 || value == 1)))
        error(['cfg.constraints.enable_xyz_snap_limit must be a scalar ' ...
            'logical value.']);
    end
    enabled = logical(value);
end
end
