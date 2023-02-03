

% 20171101 - new version of the Hilbert mask to recover "analytic" signal
% from a 3D matrix, for use in the Stockwell Transform. NPH.

%% Background
% In essence, all we are doing is creating a mask to double
% selected frequencies and set their complex conjugate pairs to zero. We
% leave and coefficients not in a complex conjugate pair alone.
% On the surface, this is fairly simple. We set all +/-ve X,Y and only
% +ve Z freqs to be doubled, and then all +/-ve X,Y and -ve Z freqs to be
% zero. However, for a fftshifted 3D fourier spectrum, there are
% 4 interlocking 2D sheets which need special attention. Where 3 planes
% intersect at one location, contain the "zeroth" frequency points, 
% ie ones not in a complex conjugate pair. For odd and even length
% dimensions, these take funny shapes.

% For an FFTSHIFTED spectrum, these are the planes that are a little bit
% special.
% 
% For a spectrum with [EVEN EVEN EVEN] dimensions:
%   .    .
%   |\   |\
%   | .----.----.           % Front, left side, base
%   .-|\-.-|\-. |           % Mid-horz plane (always a feature)
%   |\| .----.----.
%   | x-|--x-|--. |
%   .-|\|-.|\|  |\|
%    \| x----x----.
%     x-|--x-|--. |
%      \|   \|   \|
%       x----x----.
% 
% For a spectrum with [ODD ODD ODD] dimensions:
%   .    .
%   |\   |\
%   | .----.----.           
%   .-|\-.-|\-. |           % Mid-horz plane (always a feature)
%   |\| .----.----.
%   | .-|--x-|--. |
%   .-|\|-.|\|  |\|
%    \| .----.----.
%     .-|--.-|--. |
%      \|   \|   \|
%       .----.----.
% 
% "x" marks the possible zeroth frequencies (not in conj pair)
%
% Fortunately, each of these planes can follow the same pattern of masking
% depending on its odd/even dimensions. Interestingly, each plane follows
% the same 2D masking pattern as would be used for the 2DST. From what I
% can gather, when the fft encounters (for ND > 1) even-numbered
% dimensions, it dumps some complex conjugate pairs in a line at the edge
% of the spectrum. For only odd dimensions, the each coeff sits happily with
% its complex conjugate pair opposite it as a reflection through the very
% centre of the spectrum, where the 0th freq component sits. If any
% dimension is even, the fft dumps the extra coeffs in a line at the edge,
% who sit opposite their complex conjugate pairs via a reflection with a
% 0th freq component in the centre of that extra line (see below).
% 
% The number of 2s should always equal the number of 0s.
% 
% For a dimension with ODD number of elements, simply take the mask below
% from inside the line, and for EVEN elements include the elements outside
% of the lines. You essentially seem to trim the below to suit your needs
% for each of these planes.
% 
% MASK EXAMPLE:
% EVEN n_elements case | ODD n_elements case
% 
%         -ve X     +ve X
%   2 | 0 0 0 0 2 2 2 2 2
%   2 | 0 0 0 0 2 2 2 2 2 +ve Y
%   2 | 0 0 0 0 2 2 2 2 2
%   2 | 0 0 0 0 2 2 2 2 2 
%   1 | 0 0 0 0 1 2 2 2 2
%   0 | 0 0 0 0 0 2 2 2 2
%   0 | 0 0 0 0 0 2 2 2 2
%   0 | 0 0 0 0 0 2 2 2 2 -ve Y
%   0 | 0 0 0 0 0 2 2 2 2
%   --+------------------ 
%   1 | 0 0 0 0 1 2 2 2 2 
% 
% 
% Once you've set these 4 planes correctly, the rest is fairly
% straightforward - simply set everything else above the middle horizontal
% plane to be 2 (for +ve Z freqs) and set everything below (but not the
% base) to be 0 (not -ve Zfreqs). You can of course fiddle this method to
% consider different freq combinations, but tbh it's probably easier just
% to permute you matrix to what you want and put it through this code :)
% 

%% Now the function itself!
% 
% INPUTS: sz - size() of the desired mask, either 1, 2 or 3D.
% 
% OUTPUTS: m - the mask itself, same size as input dimensions
% 
% Note, we expect an fftshifted vector/matrix in each case.
% 

function m = nph_hilbertmask(sz)

m = nan(sz); % ensure output is same size as specified

sz(sz == 1) = []; % cope with 1x8, 8x1x5 etc.

mid = fix(sz/2)+1; % find midpoint(s).

switch length(sz)
    
%   1D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1 % 1D (straight from bob stockwell)
        
        if isodd(sz)
            m(1:mid-1)      = 0;
            m(mid)          = 1;
            m(mid+1:end)    = 2;
        else
            m(1)            = 1;
            m(2:mid-1)      = 0;
            m(mid)          = 1;
            m(mid+1:end)    = 2;
        end        
%         m = [ones(1-rem(sz,2),1); 2*ones(fix((sz-1)/2),1); 1; zeros(fix((sz-1)/2),1)];
    
%   2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2 % 2D (NPH and NDS method)
        % make the 2D mask featured in the preamble:
        m = mask2D(sz);
    
%   3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3 % 3D (NPH method)
        
        % First, select all +/-ve X,Y and +ve Z freqs:
        
        m(:,:,mid(3)+1:end) = 2;
        m(:,:,1:mid(3)-1)   = 0;
        
        % Do the middle slice (always done regardless of odd/even
        % dimensions)
        m(:,:,mid(3)) = mask2D(sz([1 2]));
        
        % Now, determine what extra sides you need depending on the
        % odd/even dimensions:
        
        % Base
        if iseven(sz(3))
        m(:,:,1) = mask2D(sz([1 2]));
        end
        
        % left hand side
        if iseven(sz(2))
        m(:,1,:) = mask2D(sz([1 3]));
        end
        
        % front
        if iseven(sz(1))
        m(1,:,:) = mask2D(sz([2 3]));
        end
        
end

% Sanity check:
if numel(m == 2) ~= numel(m == 0)
    disp('Warning: Problem with the Hilbert mask.')
    return
end


end


%==========================================================================
function m2 = mask2D(dims)

mid = fix(dims/2)+1;

% the order matters!
m2 = nan(dims);

m2(:,mid(2):end)            = 2; % double half
m2(1:mid(1),1:mid(2))       = 0; % zero bottom quarter
m2(mid(1):end,1:mid(2)-1)   = 0; % zero top quarter

m2(mid(1),mid(2))           = 1; % middle 0th freq

% If both sides are odd, the 2D mask is finished now! yayyy!

% Bottom and left hand side edges?
switch double(iseven(dims(1))) / double(iseven(dims(2)))
    case Inf % [1 0] odd, even
        m2(1,:) = [zeros(1,mid(2)-1) 1 2*ones(1,mid(2)-1)];
    case 0   % [0 1] even, odd
        m2(:,1) = [zeros(1,mid(1)-1) 1 2*ones(1,mid(1)-1)]';
    case 1   % [1 1] even, even
        m2(1,:) = [1 zeros(1,mid(2)-2) 1 2*ones(1,mid(2)-2)];
        m2(:,1) = [1 zeros(1,mid(1)-2) 1 2*ones(1,mid(1)-2)]';
end

% Sanity check:
if numel(m2 == 2) ~= numel(m2 == 0)
    disp('Warning: Problem with the Hilbert mask.')
    return
end

% flipud(m2 - nph_hilb2(dims))
% figure; imagesc(1:dims(1),1:dims(2),m2);


end
%==========================================================================





%==========================================================================
% ISEVEN and ISODD functions
function logikal = iseven(x)
    logikal = mod(x,2)==0;
end
function logikal = isodd(x)
    logikal = mod(x,2)==1;
end
%==========================================================================




















