function hits = findFieldsRecursive(S, patterns)
hits = strings(0,1);
walk(S, "data");
    function walk(x, path)
        if isstruct(x)
            fn = fieldnames(x);
            for i = 1:numel(x)
                for k = 1:numel(fn)
                    p2 = path + "(" + i + ")." + fn{k};
                    walk(x(i).(fn{k}), p2);
                end
            end
        elseif iscell(x)
            for i = 1:numel(x)
                walk(x{i}, path + "{" + i + "}");
            end
        else
            % If the current path name matches patterns, record it
            lowerPath = lower(path);
            for pat = patterns
                if contains(lowerPath, pat)
                    hits(end+1,1) = path; %#ok<AGROW>
                    break
                end
            end
        end
    end
end
