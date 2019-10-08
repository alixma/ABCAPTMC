function zlevs=get_levels(S, nlevels)
zmin = min(S(:));
zmax = max(S(:));
zinc = (zmax - zmin) / nlevels;
zlevs = linspace(zmin+zinc, zmax, nlevels);
end
