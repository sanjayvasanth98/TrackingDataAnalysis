function policyTag = build_case_analysis_policy_tag(basePolicyTag, convergenceFrameStep, varargin)
%BUILD_CASE_ANALYSIS_POLICY_TAG Combine analysis settings into one cache key.

parts = strings(0, 1);
parts(end+1, 1) = "caseCheckpointVersion=3";
parts(end+1, 1) = "base=" + string(basePolicyTag);
parts(end+1, 1) = "convergenceFrameStep=" + string(convergenceFrameStep);

if mod(numel(varargin), 2) ~= 0
    error('build_case_analysis_policy_tag:NameValuePairs', ...
        'Expected name/value pairs after convergenceFrameStep.');
end

for i = 1:2:numel(varargin)
    name = string(varargin{i});
    value = varargin{i + 1};
    parts(end+1, 1) = name + "=" + string(cache_value_signature(value)); %#ok<AGROW>
end

policyTag = strjoin(parts, "|");
end
