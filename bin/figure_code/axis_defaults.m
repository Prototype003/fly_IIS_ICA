function [ handle ] = axis_defaults( handle )
%AXIS_DEFAULTS
% Applies standard settings to an axes handle
%
% Inputs:
%   handle: handle to the axes
%
% Outputs:
%   handle: handle to the axes

fontsize = 11; % FontSize
titlefontsizemultiplier = 1; % TitleFontSizeMultiplier
labelfontsizemultiplier = 1;
fontname = 'arial';

set(handle,...
    'FontSize', fontsize,...
    'TitleFontSizeMultiplier', titlefontsizemultiplier,...
    'LabelFontSizeMultiplier', labelfontsizemultiplier,...
    'FontName', fontname...
    );


end

