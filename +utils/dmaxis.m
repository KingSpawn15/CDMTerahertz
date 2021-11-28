function [vectors , units] = dmaxis( tags )
% function vector = dmaxis(tags)
% Return the vectors of each axis, like X and Y, from Digital Micrograph 
% tag structure from dmread() by Andreas Korinek.
% 
% Args: 
%    - tags (structures) GMS Tags returned 
% Outputs:
%    - vectors (cell array) contains the vectors of each Axis
%    - units (cell array) units of each Axis
%
% Example 1:
%    vectors = dmaxis(tags);
%    [X, Y] = vectors{:};
%
% Author: Kangpeng Wang, Technion, 2019, wangkangpeng@msn.com
%  

% Rename to short the length of this code
ImageData = tags.ImageList.Unnamed1.ImageData;

n = 0;
while isfield(ImageData.Dimensions,['Unnamed' num2str(n)])
    % Get the length
    size = double(getfield(ImageData.Dimensions,...
           ['Unnamed' num2str(n)],'Value'));   
    % Get the Origin, Scale of each axis
    origin = double(getfield(ImageData.Calibrations.Dimension,...
        ['Unnamed' num2str(n)],'Origin','Value'));
    scale = double(getfield(ImageData.Calibrations.Dimension,...
        ['Unnamed' num2str(n)],'Scale','Value'));
    
    % personally I don't understand why Gatan I have this weired origin
    % definition, but this only works with minus operator.
    teminate = -origin + size-1; 
    
    % Generate vector for each axis
    vectors{n+1} = (-origin.*scale):scale:(teminate.*scale); %#ok<AGROW,BDSCA>
    
    % Get Units 
    units{n+1} = char(getfield(ImageData.Calibrations.Dimension,...
        ['Unnamed' num2str(n)],'Units','Value'))';  %#ok<AGROW>
    
    n = n+1;
end

end

