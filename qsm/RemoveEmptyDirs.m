function RemoveEmptyDirs(dirTop)
%REMOVEEMPTYDIRS Removes all empty directories within the provided path, including the provided path

    items = dir(dirTop);

    % Recurse into sub directories
    for d = items'
        path = strcat(d.folder, "/", d.name);
        if d.isdir && ~(endsWith(path, '/.') || endsWith(path, '/..'))
            RemoveEmptyDirs(path)
        end
    end

    % Delete this is if it is empty
    items = dir(dirTop);
    if isempty(items)
        rmdir(dirTop)
    end

end

