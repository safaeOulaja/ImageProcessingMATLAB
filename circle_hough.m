function [h, margin] = circle_hough(b, rrange, varargin)
opts = {'same' 'normalise'};
narginchk(2, 2+length(opts));
validateattributes(rrange, {'double'}, {'real' 'positive' 'vector'});
if ~all(ismember(varargin, opts))
    error('Unrecognised option');
end

% get indices of non-zero features of b
[featR, featC] = find(b);

% set up accumulator array - with a margin to avoid need for bounds checking
[nr, nc] = size(b);
nradii = length(rrange);
margin = ceil(max(rrange));
nrh = nr + 2*margin;        % increase size of accumulator
nch = nc + 2*margin;
h = zeros(nrh*nch*nradii, 1, 'uint32');  % 1-D for now, uint32 a touch faster

% get templates for circles at all radii - these specify accumulator
% elements to increment relative to a given feature
tempR = []; tempC = []; tempRad = [];
for i = 1:nradii
    [tR, tC] = circlepoints(rrange(i));
    tempR = [tempR tR]; %#ok<*AGROW>
    tempC = [tempC tC];
    tempRad = [tempRad repmat(i, 1, length(tR))];
end

% Convert offsets into linear indices into h - this is similar to sub2ind.
% Take care to avoid negative elements in either of these so can use
% uint32, which speeds up processing by a factor of more than 3 (in version
% 7.5.0.342)!
tempInd = uint32( tempR+margin + nrh*(tempC+margin) + nrh*nch*(tempRad-1) );
featInd = uint32( featR' + nrh*(featC-1)' );

% Loop over features
for f = featInd
    % shift template to be centred on this feature
    incI = tempInd + f;
    % and update the accumulator
    h(incI) = h(incI) + 1;
end

% Reshape h, convert to double, and apply options
h = reshape(double(h), nrh, nch, nradii);

if ismember('same', varargin)
    h = h(1+margin:end-margin, 1+margin:end-margin, :);
    margin = 0;
end

if ismember('normalise', varargin)
    h = bsxfun(@rdivide, h, reshape(rrange, 1, 1, length(rrange)));
end

end

