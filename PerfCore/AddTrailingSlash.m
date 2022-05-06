function path = AddTrailingSlash(path)

if (path(end) ~= '/') & (path(end) ~= '\')
    path = [path '/'];
end