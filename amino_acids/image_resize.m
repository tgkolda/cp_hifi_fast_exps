function [width,height] = image_resize(varargin)
%IMAGE_RESIZE Resize the current figure
%
%   IMAGE_RESIZE resizes the current figure to 11" x 5", which is a good 
%   size for powerpoint slides.
%
%   IMAGE_RESIZE(W,H) resizes the current figure to W in. by H in.
%
%   [W,H] = IMAGE_RESIZE(...) also returns the width and height in inches.

% Resize (keeping position of upper left corner the same)
if nargin == 0
	width = 11;     % Width in inches
	height = 5;    % Height in inches
elseif nargin == 2
	width = varargin{1};     % Width in inches
	height = varargin{2};    % Height in inches
else
	error('Invalid number of arguments. Usage: image_resize() or image_resize(W,H)')
end
pos = get(gcf, 'Position');
delta_y = max(0, height*100 - pos(4));
set(gcf, 'Position', [pos(1) (pos(2)-delta_y) width*100, height*100]); 
