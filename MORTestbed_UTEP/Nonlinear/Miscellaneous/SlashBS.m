function [sbs] = SlashBS()

if isunix
    sbs = '/';
elseif ispc
    sbs = '\';
end

end