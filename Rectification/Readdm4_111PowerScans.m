function [DATAAll,Energy,Units] = Readdm4_111PowerScans(Path,RepNumStr)

%% Create list of files in subdirectory belonging to current repetition
FolderListing = dir([Path,'*.dm4']);

FileNum = numel(FolderListing);
for i = 1:FileNum
    BoolList(i) = not(isempty(strfind(FolderListing(i).name,['_00',RepNumStr,'_']))) && length(FolderListing(i).name) == 27;
end
FolderListing = FolderListing(BoolList);

%% Read dm4 files
for FileInd = 1:length(FolderListing)
    [tags, DATA] = dmread([FolderListing(FileInd).folder,'\', FolderListing(FileInd).name]);
    DATAAll(:,:,FileInd) = double(DATA);
    [Dimensions,Units] = dmaxis(tags);
end

Energy = Dimensions{1};

end

