function barwitherr(errors,varargin)

% Check how the function has been called based on requirements for "bar"
if nargin < 3
    % This is the same as calling bar(y)
    values = varargin{1};
    xOrder = 1:size(values,1);
else
    % This means extra parameters have been specified
    if isscalar(varargin{2}) || ischar(varargin{2})
        % It is a width / property so the y values are still varargin{1}
        values = varargin{1};
        xOrder = 1:size(values,1);
    else
        % x-values have been specified so the y values are varargin{2}
        % If x-values have been specified, they could be in a random order,
        % get their indices in ascending order for use with the bar
        % locations which will be in ascending order:
        values = varargin{2};
        [tmp xOrder] = sort(varargin{1});
    end
end

% If an extra dimension is supplied for the errors then they are
% assymetric split out into upper and lower:
if ndims(errors) == ndims(values)+1
    lowerErrors = errors(:,:,1);
    upperErrors = errors(:,:,2);
elseif isvector(values)~=isvector(errors)
    lowerErrors = errors(:,1);
    upperErrors = errors(:,2);
else
    lowerErrors = errors;
    upperErrors = errors;
end
    
% Check that the size of "errors" corresponsds to the size of the y-values.
% Arbitrarily using lower errors as indicative.
if any(size(values) ~= size(lowerErrors))
    error('The values and errors have to be the same length')
end

[nRows nCols] = size(values);
handles.bar = bar(varargin{:}); % standard implementation of bar fn
hold on

if nRows > 1
    for col = 1:nCols
        % Extract the x location data needed for the errorbar plots:
        x = get(get(handles.bar(col),'children'),'xdata');
        % Use the mean x values to call the standard errorbar fn; the
        % errorbars will now be centred on each bar; these are in ascending
        % order so use xOrder to ensure y values and errors are too:
        errorbar(mean(x,1),values(xOrder,col),lowerErrors(xOrder,col), upperErrors(xOrder, col), '.k')
    end
else
    x = get(get(handles.bar,'children'),'xdata');
    errorbar(mean(x,1),values,errors,'.k')
end

hold off