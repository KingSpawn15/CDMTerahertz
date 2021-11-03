function [DATAAll,Energy,Units] = Readdm4_Pol(Path)

%% Read dm4 files
[tags, DATA] = dmread(Path);
DATAAll = double(DATA);
[Dimensions,Units] = dmaxis(tags);
Energy = Dimensions{1};

end

