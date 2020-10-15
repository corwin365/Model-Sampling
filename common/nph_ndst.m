
% nph_ndst.m
% A 1, 2 and 3-D application of the Stockwell Transform, developed by
% Neil Hindley, University of Bath, 2017.

% CHANGE LOG ==============================================================

% - August 2018: I did it. I finally did it. Found a way to cope with the
% amplitude underestimation for higher dimensions. We now take the Hilbert
% transform of the input N-D image that has been pseudo-bandpassed by all
% the gaussian windows that we consider (see below) in the ST. As a result,
% the Hilbert Boosting method has been removed.

% - 2018 Several changes since I last updated, mostly to make it more N-D
% friendly. First, the hilbert masking step is now fully N-D, we simply
% list all the complex conjugate pairs and set one to zero double the other
% and put back. It's essentially abitrary where they were, so long as the
% pair has been addressed correctly. This meant that the gaussian windows
% had to be able to cope with arbitrary masking, so I made double
% gaussians, where each has a reflection in the opposite quadrant across
% the zero frequency line. I then stopped any leakage into adjacent
% quadrants, as sometimes I was getting a better signal from a quadrant I
% wasn't windowing when collapsing the spectrum. I've also swapped the
% repmatting in the main loop for ndgrids. Oh and I now assemble the
% gaussians outside the main loop as vectors and then grid, mulitply and
% add them accordingly within the loop, minimising cost. Also this step is now
% full N-D, it's just really the main loop now that has switch cases. Oh
% also the decision to make the masking step abitrary was so that it didn't
% matter what positive/negative scales combinations the user put in. This
% is much, much more robust now. 
% EDIT: I've held off on the spillage thing for now. Maybe I'll uncomment
% it back in at some point.

% - late 2017 new method for computing gaussian windows takes a fraction of the time.
% method works by taking each dimension's exponential as a vector then
% repmatting and multiplying in the loop, thereby removing the exponential.
%
% - 20171218 "Hilbert Boosting": Use a method that convolves the
% reconstructed ST spectrum with the Hilbert Transform of the original data
% to get a better estimate of what the original amplitude was. This helps
% take into account the packet-like nature of waves.
%
% - 20180705 New N-D method for computing the Hilbert mask, and "Double
% Gaussian Windows" are now the norm. OMG this method is simpler!!!!



%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ST = nph_ndst(IN,scales,point_spacing,c,varargin)

% *nomenclature: N = number of dimensions, n = number of elements in a
% particular dimension.
%
% IN [REQUIRED] - N-D vector/matrix to be transformed. 1, 2, and 3-D supported.
%
% scales - 1xN cell object, where N is dimensions, each containing a 1xm
% length vector of integer scales with which to analyse the relevant
% dimenion. These correspond to integer fractions of the maximum possible
% wavelength, ie the full length of the input series, up to the Nyquist.
% For 1-D, just entering a vector is fine, no need for a cell.
% Scales should be anything withing the range [-(n/2-1):-1 1:1:(n/2-1)],
% where n is the number of elements in a particular dimension. Defaults are
% roughly [-(n/3):(n/3)] for each dimension for the 2DST and 3DST, and
% 1:Nyquist (1:n-1) for the 1DST.
%
% point_spacing - 1xN vector, containing the real physical separation
% between points, be it time or distance. Must have regular sampling
% intervals.
%
% c - scaling parameter for the Guassian windows. See Hindley et al., AMT
% (2016) for details on the effect of this.
%
% [OPTIONAL ARGUMENTS]
%
% 'quick' (default) | 'full' - choose whether to output the full
% 2N dimensional S-transform complex output. For 1-D, default is 'full', as
% it's only small, but for 2- and 3-D, default is 'quick', which provides
% the dominant spectral component at each location in the data. Remember,
% the 3-D output is 6-D, which can be pretty massive, so use with caution.
%
% 'zeromean' (default) | 'nozeromean' - choose whether to compute the
% NDST on the zero-mean signal or not. Should really be doing this as a
% rule anyway.
%
% % % % % REMOVED: 'boost', new hilbert boosting is done by default.
% % % % % 'boost' (default for 2D/3D) | 'noboost' (default for 1D) - choose
% % % % % whether to apply the Hilbert Boosting method described below to try to
% % % % % get a better estimate of wave-packet amplitudes, which are typically
% % % % % underestimated by the 2D and 3D ST, but less so by the 1D ST.
%
% [OUTPUTS]
%
% ST.ST - (m{1,2,3} x n{1,2,3}) Stockwell Tranform complex cospectrum.
%
% ST.C - N-dimensional complex cospectrum of the dominant (largest spectral
% amplitude) frequency at each location. This will either be boosted/not
% boosted depending on whether you want to boost the amplitude to cope with
% the packet-like nature of waves.
%
% ST.A - abs() of ST.C, the instantaneous amplitude of the dominant
% frequency at each location.
%
% ST.R - real() part of ST.C, provides "reconstruction" of the wavefield as
% the NDST saw it.
%
% ST.F{1,2,3} - dominant N-dimensional spatial frequencies at each location
% (inverse of distance or time, no 2*pi)
%
% ST.freqs - 1xN cell object of the frequencies that were analysed,
% computed from the scales that were inputted.
%
% ST.HA - Absolute part of the Hilbert Transform of the psuedo-bandpassed
% data. All the Gaussian windows are combined and applied to the FFT
% spectrum, justa  rough bandpassed filter, then the Hilbert transform is
% applied to this to find the phase-invariant amplitude at each location.
% This is superior because it is the amplitude of only those wavelengths
% that we have considered in the S-transform.
% 
% ST.HR - Real part of the above.
% 
% 
% % % % % % REMOVED:
% % % % % % ST.BoostFactor - factor by which ST.C has been boosted by the Hilbert
% % % % % % Boosting method.
%
%


function ST = nph_ndst(IN,varargin)

% % figure; hold all;

IN = squeeze(IN);

%% Determine if 1D 2D 3D input ============================================
sz = size(IN);
sz = sz(sz ~= 1);
type = length(sz);

% fix the annoying 1xN versus Nx1 problem:
if type == 1
    IN = reshape(IN,[1 length(IN)]);
end
osz = size(IN); % original corrected size in

%% VARARGIN ===============================================================
% NEW PLAN: Split varargin into inputs and string options:
inputs = {};
options = {};
for v = 1:length(varargin)
    if ischar(varargin{v}) || isstring(varargin{v})  % options must be a character string
        options{v} = varargin{v};
    else
        inputs{v} = varargin{v};
    end
end

% Default SCALES:
default_scales = cell(1,type);
for i = 1:type
    if i < type
        default_scales{i} = -15:15;
    else
        default_scales{i} = 1:15;
        % just positive freqs for the last dimension, save repeating yourself.
    end
end

% Default POINT SPACING:
default_point_spacing = ones(1,type);

% Default scaling parameter C:
default_c = 0.25 * ones(1,type);

% Assign defaults if they're missing:
switch length(inputs)
    case 0 % no inputs
        scales = default_scales;
        point_spacing = default_point_spacing;
        c = default_c;
    case 1 % just scales
        scales = inputs{1};
        point_spacing = default_point_spacing;
        c = default_c;
    case 2 % just scales and point_spacing:
        scales = inputs{1};
        point_spacing = inputs{2};
        c = default_c;
    otherwise % scales, point_spacing and c:
        scales = inputs{1};
        point_spacing = inputs{2};
        c = inputs{3};
end

% Specify some defaults for options:
switch type
    case 1
        fullflag = 1;
        boostflag = 0;
        zeromeanflag = 1;
    otherwise
        fullflag = 0;
        boostflag = 1;
        zeromeanflag = 1;
end

% Now overwrite these if the user requests it:

% FULL 2-D, 4-D, 6-D etc S-TRANSFORM OBJECT
% Determine if we want full 2D, 4D, 6D etc output:
if any(strcmpi(options,'full'))
    fullflag = 1;
    % warning for 3DST:
    if type == 3
        if ~isempty(inputs) % if scales are inputted
            s = inputs{1};
            w = whos('IN');
            Mb = (w.bytes ./ 1024 ./ 1024);
            switch class(IN)
                case 'double'
                    Mb = Mb ./ 2;
            end
            Mb = Mb .* length(s{1})*length(s{2})*length(s{3});
            warning(['Outputing a full 6-D S-transform object at single precision. This might require up to ' num2str(Mb) ' Mb of memory.'])
        end
    end
end

% ZERO MEAN
% Determine if we want to use zero-mean (as a rule we should really get
% in the habit of giving signals zero-mean before the NDST. This means
% of course that any underlying trend will manifest as a signal, but it
% results in fairer recovery of stuff which sits upon that trend.
if any(strcmpi(options,'nozeromean'))
    zeromeanflag = 0;
end
if any(strcmpi(options,'zeromean'))
    zeromeanflag = 1;
end

% HILBERT BOOSTING
% Determine if we want to use the "Hilbert Boosting" method to try to
% more accurately measure the amplitudes of wave packets, which are
% always underestimated by the ST.
if any(strcmpi(options,'noboost'))
    boostflag = 0;
end
if any(strcmpi(options,'boost'))
    boostflag = 1;
end


%% PARSE INPUT SIZES ======================================================

switch type
    case 1
    otherwise
        if ~iscell(scales) || length(scales) ~= type
            error(['Error: Scales must be 1x' num2str(type) ' cell of 1D vectors for each dimension.'])
        end
        if ~isnumeric(point_spacing) || length(point_spacing) ~= type
            error(['Error: Point spacing must be a 1x' num2str(type) ' vector.'])
        end
        if ~isnumeric(c) || length(c) ~= type
            error(['Error: Scaling parameter "c" must be a 1x' num2str(type) ' vector.'])
        end
end

% for 1D, put scales as a cell just for ease of coding later
if type == 1 && ~iscell(scales), scales = {scales}; end

%% PARSE FOR NON-INTEGER SCALES ===========================================

for i = 1:length(scales) % integer scales, not equal to zero.
    sc = scales{i};
    sc = unique(fix(sc),'stable'); % unique sorts it by default :(
    scales(i) = {sc(sc ~= 0)};
end

%% INITIALISE OUTPUTS =====================================================

% initialise:
ST = struct;

% record inputs:
ST.IN = IN;
ST.scales = scales;
ST.point_spacing = point_spacing;
ST.c = c;
ST.AmplitudeBoosting = boostflag;

if ~isempty(options), ST.Options = options; end

% sort out scales and frequencies
numscales = zeros(1,type);
for t = 1:type
    numscales(t) = length(scales{t});
    ST.freqs{t} = scales{t} ./ (sz(t)*point_spacing(t));
end

% for 1D, un-cell scales and freqs for ease of reading in the outputs
if type == 1 && iscell(ST.freqs)
    ST.scales = ST.scales{1};
    ST.freqs = ST.freqs{1};
end

% initialise main outputs:

switch type
    case 1 % 1DST always has the full spectrum
        ST.ST = zeros([numscales sz],'single');
        ST.C = zeros(osz,'single');
        ST.F1 = zeros(osz,'single');
    case 2
        if fullflag
            ST.ST = zeros([numscales(1) numscales(2) sz],'single');
            ST.C = zeros(osz,'single');
        else
            ST.C = zeros(osz,'single');
        end
        
        ST.F1 = zeros(osz,'single');
        ST.F2 = zeros(osz,'single');
        
    case 3
        
        if fullflag
            ST.ST = zeros([numscales sz],'single');
            ST.C = zeros(osz,'single');
        else
            ST.C = zeros(osz,'single');
        end
        
        ST.F1 = zeros(osz,'single');
        ST.F2 = zeros(osz,'single');
        ST.F3 = zeros(osz,'single');
        
end

%% ZERO-MEAN ==============================================================
if zeromeanflag
    IN = IN - mean(IN(:));
end

%% STEP 1: FFTN ===========================================================

F = single(fftn(IN));

% dc_comps = imag(F) == 0;

%% STEP 2: GATHER FFT WISDOM ==============================================
% gather wisdom for the IFFTN() in the loop below...
if type > 1 % don't bother for the 1DST
    ST.fftwisdom = gatherfftwisdom(F,'ifftn');
    fftw('swisdom',ST.fftwisdom); % Apply the wisdom for SINGLE precision.
end

%% STEP 3: HILBERT MASK ===================================================

FM = F .* nph_mask(F);

%% STEP 4: ASSEMBLE GAUSSIAN WINDOWS ======================================
% now fully ND compatible, pre-make gaussian vectors into a structure

GW = struct;

for i = 1:length(scales)
    
    % Define Gaussian vector storage structure:
    GW(i).gvecA = nan(length(scales{i}),sz(i));
    GW(i).gvecB = nan(length(scales{i}),sz(i));
    
    for j = 1:length(scales{i})
        
        % First, make a coord system in fftshifted space (easier to work in)
        switch iseven(sz(i))
            case 1
                N = (sz(i)/2)-1;
                x = -(N+1):N;
%                 x = [0 -N:N];
            case 0
                N = (sz(i)-1)/2;
                x = -N:N;
        end
        
        % ifftshift:
        x = ifftshift(x);
        
        % Evaluate Gaussians in the normal way:
        % normal window:
%         gwA = exp( (-2*pi^2) * (c(i)/scales{i}(j))^2 * ((x.^2) - 2*x*scales{i}(j) + scales{i}(j)^2));
        gwA = exp( (-2*pi^2) * (c(i)/scales{i}(j))^2 * (x - scales{i}(j)).^2 );
        
        % a mirror image window:
%         gwB = exp( (-2*pi^2) * (c(i)/scales{i}(j))^2 * ((x.^2) + 2*x*scales{i}(j) + scales{i}(j)^2));
        gwB = exp( (-2*pi^2) * (c(i)/scales{i}(j))^2 * (x + scales{i}(j)).^2 );
        
%         % Remove any spillage into adjacent quadrants:
%         switch sign(scales{i}(j))
%             case 1
%                 gwA(x < 0) = 0;
%                 gwB(x > 0) = 0;
%             case -1
%                 gwA(x > 0) = 0;
%                 gwB(x < 0) = 0;
%         end
        
        % and assign:
        GW(i).gvecA(j,:) = gwA;
        GW(i).gvecB(j,:) = gwB;
    end
    
end

%% MAKE A GAUSSIAN WINDOW STORE to find what freqs we computed:
gws = zeros(osz);

%% STEP 5: APPLY THE GAUSSIAN WINDOWS AND IFFTN ===========================
% figure; hold all; figpos([1 0.5])
% rows = 1; cols = 5;
% subplot(rows,cols,1);
% imagesc(ST.IN); ydir; axis tight;
% subplot(rows,cols,2);
% imagesc(abs(fftshift(F))); ydir; axis tight;

switch type
    case 1 % 1DST
%         for i1 = 1:length(scales{1})
        for i1 = randperm(length(scales{1}))    
            % for this gaussian window...
            gwA = GW(1).gvecA(i1,:);
            gwB = GW(1).gvecB(i1,:);
            gw = gwA + gwB;
            % gaussian window storage for hilbert amplitude later:
            gws = gws + gw;
            % Always compute the full 2D spectrum, it's not much :)
            FM_voice = ifftn(FM .* gw);
            % insert into S-transform
            ST.ST(i1,:) = FM_voice;
            % in any case, do the collapsed rapide spectrum too:
            loc = abs(FM_voice) > abs(ST.C);
            ST.C(loc) = FM_voice(loc);
            ST.F1(loc) = ST.freqs(i1);
        end
    case 2 % 2DST
%         for i1 = 1:length(scales{1})
%             for i2 = 1:length(scales{2})
        for i1 = randperm(length(scales{1}))
            for i2 = randperm(length(scales{2}))        
                % assemble gaussian window
                [gw1A,gw2A] = ndgrid(GW(1).gvecA(i1,:),GW(2).gvecA(i2,:));
                [gw1B,gw2B] = ndgrid(GW(1).gvecB(i1,:),GW(2).gvecB(i2,:));
                
                gw = (gw1A .* gw2A) + (gw1B .* gw2B);
                
%                 gw(dc_comps) = 1;
                
                % gaussian window storage for hilbert amplitude later:
                gws = gws + gw;
                
%                 subplot(rows,cols,3);
%                 cla; imagesc(fftshift(gw)); ydir;
%                 clim([0 1]); axis tight;
%                 ylabel(scales{1}(i1));
%                 xlabel(scales{2}(i2));

                % apply it
                FM_voice = ifftn(FM .* gw);
                
                if fullflag % full 4D spectrum, if required
                    ST.ST(i1,i2,:,:) = FM_voice;
                end
                % in any case, do the collapsed rapide spectrum
                loc = abs(FM_voice) > abs(ST.C);
                ST.C(loc) = FM_voice(loc);
                ST.F1(loc) = ST.freqs{1}(i1);
                ST.F2(loc) = ST.freqs{2}(i2);
                
%                 subplot(rows,cols,4);
%                 cla; imagesc(1./ST.F1); ydir;
%                 clim(minmax(1./ST.freqs{1})); axis tight;
%                 title('1/F1')
%                 
%                 subplot(rows,cols,5);
%                 cla; imagesc(1./ST.F2); ydir;
%                 clim(minmax(1./ST.freqs{2})); axis tight;
%                 title('1/F2')
%                 
%                 drawnow;
                
            end
        end
    case 3 % 3DST
%         for i1 = 1:numscales(1)
%             for i2 = 1:numscales(2)
%                 for i3 = 1:numscales(3)
        for i1 = randperm(length(scales{1}))
            for i2 = randperm(length(scales{2}))
                for i3 = randperm(length(scales{3}))
                    % assemble gaussian window
                    [gw1A,gw2A,gw3A] = ndgrid(GW(1).gvecA(i1,:),GW(2).gvecA(i2,:),GW(3).gvecA(i3,:));
                    [gw1B,gw2B,gw3B] = ndgrid(GW(1).gvecB(i1,:),GW(2).gvecB(i2,:),GW(3).gvecB(i3,:));
                    gw = (gw1A .* gw2A .* gw3A) + (gw1B .* gw2B .* gw3B);
                    % gaussian window storage for hilbert amplitude later:
                    gws = gws + gw;
                    % apply it
                    FM_voice = ifftn(FM .* gw);
                    if fullflag % full 6D spectrum, if required
                        ST.ST(i1,i2,i3,:,:,:) = FM_voice;
                    end
                    % in any case, do the collapsed rapide spectrum
                    loc = abs(FM_voice) > abs(ST.C);
                    ST.C(loc) = FM_voice(loc);
                    ST.F1(loc) = ST.freqs{1}(i1);
                    ST.F2(loc) = ST.freqs{2}(i2);
                    ST.F3(loc) = ST.freqs{3}(i3);
                end
            end
        end
end

%% GET ABS AND REAL PARTS =================================================
% for the lazy (yes you corwin)
ST.A = abs(ST.C);
ST.R = real(ST.C);


%% (ST-FILTERED) HILBERT AMPLITUDE ========================================
% Replaces the old (well not that old) Hilbert Boosting method.

% take all those Gaussian windows from earlier and normalise to 1:
gws = gws ./ max(gws(:));
gws(gws < 0) = 0; % check for anomalies.
% this gives us a nice blob showing which parts of the fft spectrum we've
% considered in this ST, given the input scales.

% Now use this like a filter on your input data, and take the Hilbert
% transform of the result:
H = ifftn(FM .* gws);

% Assign the real() and abs() amplitudes:
ST.HA = abs(H);
ST.HR = real(H);

ST.allGWs = gws;

% and you're done. Simple as that. I think back to all the years I've been
% thinking of how to do this, and here we are in a few lines. Mad.


% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % EDIT: Corwin didn't like the Hilbert amplitude, so let's try and do a
% % % % convolution again but using the "filtered" spectrum.
% % % 
% % % FR = fftn(ST.R);
% % % FH = fftn(ST.HR);
% % % 
% % % meanfactor = nanmean(abs(FH(:))) / nanmean(abs(FR(:)));
% % % 
% % % cv = sqrt(FH .* conj(FR .* meanfactor));
% % % 
% % % % cv = cv .* nph_mask(F);
% % % 
% % % disp('fart')
% % % 
% % % % meanfactor = nanmean(abs(ST.HA(:))) / nanmean(abs(ST.A(:)));
% % % % 
% % % % cv = sqrt(H .* conj(ST.C .* meanfactor));
% % % % 
% % % % cv = cv .* nanmean(abs(H(:)))/nanmean(abs(cv(:)));
% % % % 
% % % % ST.covA = abs(cv);
% % % % ST.covR = real(cv);
% % % 
% % % % % % % % Fgw = F .* gws;
% % % % 
% % % % Fst = fftn(ST.R);
% % % % 
% % % % factor = nanmean(abs(Fgw(:))) / nanmean(abs(Fst(:)));
% % % % 
% % % % cv = sqrt(Fgw .* conj(Fst .* factor));
% % % % 
% % % % C = ST.C = ST.C .* (abs(cv) ./ abs(ST.C));
% % Fw = F .* gws;
% 
% % Hf = FM .* gws;
% hmean = nanmean(abs(Fw(:)));
% stmean = nanmean(abs(ifftn(ST.R)));
% 
% C = ST.C .* (hmean/stmean);
% 
% cv = sqrt(C .* conj(Hf));
% 
% ST.covA = abs(cv);
% ST.covR = real(cv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








% % % % %
% % % % % %% STEP 4: GAUSSIAN WINDOW COORDINATE SYSTEM ==============================
% % % % % % create it outside the loop, to save processing time...
% % % % % switch type
% % % % %     % NEW GAUSSIAN WINDOW METHOD:
% % % % %     case 1
% % % % %         V1 = single((1:sz(1)) - (fix(sz(1)/2)+1));
% % % % %         V1 = ifftshift(V1); % do the ifftshift here, save it from the loop
% % % % %         V1_sq = V1.^2; % do the squared coord frame here, save doing it in the loop
% % % % %     case 2
% % % % %         V1 = single((1:sz(1)) - (fix(sz(1)/2)+1))'; % transpose is important
% % % % %         V2 = single((1:sz(2)) - (fix(sz(2)/2)+1));
% % % % %         V1 = ifftshift(V1);
% % % % %         V2 = ifftshift(V2);
% % % % %         V1_sq = V1.^2;
% % % % %         V2_sq = V2.^2;
% % % % %     case 3
% % % % %         V1 = single((1:sz(1)) - (fix(sz(1)/2)+1))'; % transpose is important
% % % % %         V2 = single((1:sz(2)) - (fix(sz(2)/2)+1));
% % % % %         V3 = single((1:sz(3)) - (fix(sz(3)/2)+1));
% % % % %         V3 = reshape(V3,[1 1 length(V3)]); % reshape for easy repmatting later.
% % % % %
% % % % %         V1 = ifftshift(V1);
% % % % %         V2 = ifftshift(V2);
% % % % %         V3 = ifftshift(V3);
% % % % %         V1_sq = V1.^2;
% % % % %         V2_sq = V2.^2;
% % % % %         V3_sq = V3.^2;
% % % % % end
% % % % %
% % % % % %% STEP 5: APPLY GAUSSIAN WINDOW AND IFFTN ================================
% % % % % % go!
% % % % % switch type
% % % % %     %==========================================================================
% % % % %     case 1 % 1DST
% % % % %
% % % % %         for i1 = 1:numscales(1), freq1 = scales{1}(i1);
% % % % %             A1 = exp( (-2*pi^2) * (c(1)/abs(freq1))^2 * (V1_sq - 2*V1*freq1 + freq1^2));
% % % % %             B1 = exp( (-2*pi^2) * (c(1)/abs(freq1))^2 * (V1_sq + 2*V1*freq1 + freq1^2));
% % % % %             gw = A1 + B1;
% % % % %
% % % % %             FM_voice = ifft(FM .* gw);
% % % % %
% % % % %             ST.ST(i1,:) = FM_voice;
% % % % %             % Always compute the full 2D spectrum, it's not much :)
% % % % %
% % % % %             % in any case, do the collapsed rapide spectrum too:
% % % % %             loc = abs(FM_voice) > abs(ST.C);
% % % % %             ST.C(loc) = FM_voice(loc);
% % % % %             ST.F1(loc) = ST.freqs(i1);
% % % % %
% % % % %         end
% % % % %
% % % % %         %         [~,inds] = max(abs(ST.ST),[],1);
% % % % %         %         ST.F1 = ST.freqs(inds);
% % % % %
% % % % %         ST.A = abs(ST.C);
% % % % %         ST.R = real(ST.C);
% % % % %
% % % % %         %==========================================================================
% % % % %     case 2 % 2DST
% % % % %         % assemble gaussian components only when you have to:
% % % % %         for i1 = 1:numscales(1), freq1 = scales{1}(i1);
% % % % %             A1 = exp( (-2*pi^2) * (c(1)/abs(freq1))^2 * (V1_sq - 2*V1*freq1 + freq1^2));
% % % % %             B1 = exp( (-2*pi^2) * (c(1)/abs(freq1))^2 * (V1_sq + 2*V1*freq1 + freq1^2));
% % % % %             for i2 = 1:numscales(2), freq2 = scales{2}(i2);
% % % % %                 A2 = exp( (-2*pi^2) * (c(2)/abs(freq2))^2 * (V2_sq - 2*V2*freq2 + freq2^2));
% % % % %                 B2 = exp( (-2*pi^2) * (c(2)/abs(freq2))^2 * (V2_sq + 2*V2*freq2 + freq2^2));
% % % % %
% % % % %                 G1 = repmat(A1,1,osz(2));
% % % % %                 G2 = repmat(A2,osz(1),1);
% % % % %                 G1b = repmat(B1,1,osz(2));
% % % % %                 G2b = repmat(B2,osz(1),1);
% % % % %                 gw = (G1 .* G2) + (G1b .* G2b);
% % % % %
% % % % %                 FM_voice = ifftn(FM .* gw);
% % % % %
% % % % %                 if fullflag % full 4D spectrum, if required
% % % % %                     ST.ST(i1,i2,:,:) = FM_voice;
% % % % %                 end
% % % % %
% % % % %                 % in any case, do the collapsed rapide spectrum
% % % % %                 loc = abs(FM_voice) > abs(ST.C);
% % % % %                 ST.C(loc) = FM_voice(loc);
% % % % %                 ST.F1(loc) = ST.freqs{1}(i1);
% % % % %                 ST.F2(loc) = ST.freqs{2}(i2);
% % % % %
% % % % %             end
% % % % %         end
% % % % %
% % % % %         ST.A = abs(ST.C);
% % % % %         ST.R = real(ST.C);
% % % % %
% % % % %         %==========================================================================
% % % % %     case 3 % 3DST
% % % % %         % assemble gaussian components only when you have to:
% % % % %         for i1 = 1:numscales(1), freq1 = scales{1}(i1);
% % % % %             A1 = exp( (-2*pi^2) * (c(1)/abs(freq1))^2 * (V1_sq - 2*V1*freq1 + freq1^2));
% % % % %             B1 = exp( (-2*pi^2) * (c(1)/abs(freq1))^2 * (V1_sq + 2*V1*freq1 + freq1^2));
% % % % %             for i2 = 1:numscales(2), freq2 = scales{2}(i2);
% % % % %                 A2 = exp( (-2*pi^2) * (c(2)/abs(freq2))^2 * (V2_sq - 2*V2*freq2 + freq2^2));
% % % % %                 B2 = exp( (-2*pi^2) * (c(2)/abs(freq2))^2 * (V2_sq + 2*V2*freq2 + freq2^2));
% % % % %                 for i3 = 1:numscales(3), freq3 = scales{3}(i3);
% % % % %                     A3 = exp( (-2*pi^2) * (c(3)/abs(freq3))^2 * (V3_sq - 2*V3*freq3 + freq3^2));
% % % % %                     B3 = exp( (-2*pi^2) * (c(3)/abs(freq3))^2 * (V3_sq + 2*V3*freq3 + freq3^2));
% % % % %
% % % % %                     G1  = repmat(A1,1,osz(2),osz(3));
% % % % %                     G2  = repmat(A2,osz(1),1,osz(3));
% % % % %                     G3  = repmat(A3,osz(1),osz(2),1);
% % % % %                     G1b = repmat(B1,1,osz(2),osz(3));
% % % % %                     G2b = repmat(B2,osz(1),1,osz(3));
% % % % %                     G3b = repmat(B3,osz(1),osz(2),1);
% % % % %                     gw = (G1 .* G2 .* G3) + (G1b .* G2b .* G3b);
% % % % %
% % % % %                     FM_voice = ifftn(FM .* gw);
% % % % %
% % % % %                     if fullflag % full 6D spectrum, if required
% % % % %                         ST.ST(i1,i2,i3,:,:,:) = FM_voice;
% % % % %                     end
% % % % %
% % % % %                     % in any case, do the collapsed rapide spectrum
% % % % %                     loc = abs(FM_voice) > abs(ST.C);
% % % % %                     ST.C(loc) = FM_voice(loc);
% % % % %                     ST.F1(loc) = ST.freqs{1}(i1);
% % % % %                     ST.F2(loc) = ST.freqs{2}(i2);
% % % % %                     ST.F3(loc) = ST.freqs{3}(i3);
% % % % %
% % % % %                 end
% % % % %             end
% % % % %         end
% % % % %
% % % % %         ST.A = abs(ST.C);
% % % % %         ST.R = real(ST.C);
% % % % %
% % % % %         %==========================================================================
% % % % %
% % % % % end

%% HILBERT BOOSTING =======================================================
if boostflag
    % %% HILBERT AMPLITUDE BOOSTING ALGORITHM:
    %%%%% NEW!!!! %%%%%
    % 2 - take hilbert transform (apply hilbert mask in fft space)
    % 3 - take complex hilbert spectrum and ST spectrum, then boost the ST
    % spectrum by the fraction of their abs() means, such that they have the
    % same abs() mean.
    % 4 - take the sqrt() of the product of the ST spectrum and the conj() of
    % the hilbert spectrum to get covarying amplitude.
    
    C_orig = ST.C;
    
    % Get complex "Hilbert" spectrum of instantaneous amplitudes:
    %     H_in = ifftn(ifftshift(fftshift(fftn(IN)) .* nph_hilbertmask(size(IN))));
    F = fftn(IN);
    H_in = ifftn(F .* nph_mask(F));
    % note - not technically hilbert transform, that's only the complex part of
    % the complex instantaneous phase object, where the real part is the
    % original signal. Look up defintion of HIlbert transform.
    
    % set both to have the same mean, so that the spectral energy is the same
    % just redistributed in different frequencies:
    hmean = nanmean(abs(H_in(:)));
    stmean = nanmean(abs(ST.C(:)));
    ST.C = ST.C .* (hmean/stmean);
    
    % Covary the two to get covarying amplitude between them:
    cv = sqrt(ST.C .* conj(H_in));
    ST.C = ST.C .* (abs(cv) ./ abs(ST.C)); % boost by fractional difference at each location
    
    ST.A = abs(ST.C);
    ST.R = real(ST.C);
    
    % NOTE this does NOT apply to the full ST object (the big 2-D/4-D/6-D
    % one). The reason for this is that, of course, the instantaneous
    % amplitude from the Hilbert transform is not defined for when the
    % signal is decomposed for each frequency voice as in the Stockwell
    % transform otherwise it would be, well, a Stockwell transform.
    
    ST.BoostFactor = abs(ST.C) ./ abs(C_orig);
    
end









%==========================================================================
% NDST FIN
%==========================================================================
end




%==========================================================================
% NESTED FUNCTIONS
%==========================================================================

% HILBERT MASK VERSION 2 ==================================================
% A newer, and much more simple, N-D approach to obtaining the analytic
% signal.
% We simply do what we say we do in the paper: find all the
% complex-conjugate pairs, set one to be zero and double the other. All
% coefficients not in a complex-conjugate pair are left unchanged.
% This approach uses linear indeces, so is N-D! Easy!

% Inputs: the fourier spectrum of the input data. Doesn't matter whether
% you've fftshift-ed it or not, the mask will be based on whatever
% arrangement you feed in.

function m = nph_mask(F)

% mask template of NaNs
m = nan(size(F));

% first, put in the zero freqs:
m(imag(F) == 0) = 1;

% now, find pairs:
[a,ib] = sort(abs(imag(F(:))));

% get rid of the zeros before reshaping:
ib = ib(a ~= 0);
a = a(a ~= 0);

% now reshape:
% this should always work - after the zeros are taken out there should
% always be an even number of complex conjugate pairs remaining.
% note: you need the transpose ' here due to the way reshape re-lists things.
ar = reshape(a',2,length(a)/2);
ibr = reshape(ib',2,length(ib)/2);

% now assign 2s and 0s (doesn't matter which order):
m(ibr(1,:)) = 2;
m(ibr(2,:)) = 0;

% and you're done!!

end



% % % % % % OLD HILBERT MASK (now disused, but useful background reading) ===========
% % % % % %
% % % % % % Matlab's hilb.m function only computes the transform across rows for N-D
% % % % % % matirices, need to write our own:
% % % % % %
% % % % % % m = nph_hilbertmask(sz,varargin);
% % % % % %
% % % % % % CHANGE LOG:
% % % % % %
% % % % % % 20171101 - new version of the Hilbert mask to recover "analytic" signal
% % % % % % from a 3D matrix, for use in the Stockwell Transform. NPH.
% % % % % %
% % % % % % 20180311 - Added an input 'neg' flag for the 3DST to select only negative
% % % % % % z freqs. We could have computed this based on inputting the individual
% % % % % % scale, but this would have meant computing the mask every timestep.
% % % % % %
% % % % % %
% % % % % %%% Background
% % % % % % In essence, all we are doing is creating a mask to double
% % % % % % selected frequencies and set their complex conjugate pairs to zero. We
% % % % % % leave and coefficients not in a complex conjugate pair alone.
% % % % % % On the surface, this is fairly simple. We set all +/-ve X,Y and only
% % % % % % +ve Z freqs to be doubled, and then all +/-ve X,Y and -ve Z freqs to be
% % % % % % zero. However, for a fftshifted 3D fourier spectrum, there are
% % % % % % 4 interlocking 2D sheets which need special attention. Where 3 planes
% % % % % % intersect at one location, they contain the "zeroth" frequency points,
% % % % % % ie ones not in a complex conjugate pair. For odd and even length
% % % % % % dimensions, these take funny shapes. This is just a quirk of the
% % % % % % arrangement of the frequencies in the FFT spectrum to result in an
% % % % % % output that is the same size as the input.
% % % % %
% % % % % % For an FFTSHIFTED spectrum, these are the 4 planes that are a little bit
% % % % % % special.
% % % % % %
% % % % % % For a spectrum with [EVEN EVEN EVEN] dimensions:
% % % % % %   .    .
% % % % % %   |\   |\
% % % % % %   | .----.----.           % Front, left side, base
% % % % % %   .-|\-.-|\-. |           % Mid-horz plane (always a feature)
% % % % % %   |\| .----.----.
% % % % % %   | x-|--x-|--. |
% % % % % %   .-|\|-.|\|  |\|
% % % % % %    \| x----x----.
% % % % % %     x-|--x-|--. |
% % % % % %      \|   \|   \|
% % % % % %       x----x----.
% % % % % %
% % % % % % For a spectrum with [ODD ODD ODD] dimensions:
% % % % % %   .    .
% % % % % %   |\   |\
% % % % % %   | .----.----.
% % % % % %   .-|\-.-|\-. |           % Mid-horz plane (always a feature)
% % % % % %   |\| .----.----.
% % % % % %   | .-|--x-|--. |
% % % % % %   .-|\|-.|\|  |\|
% % % % % %    \| .----.----.
% % % % % %     .-|--.-|--. |
% % % % % %      \|   \|   \|
% % % % % %       .----.----.
% % % % % %
% % % % % % "x" marks the possible zeroth frequencies (not in conj pair)
% % % % % %
% % % % % % Fortunately, each of these planes can follow the same pattern of masking
% % % % % % depending on its odd/even dimensions. Interestingly, each plane follows
% % % % % % the same 2D masking pattern as would be used for the 2DST. From what I
% % % % % % can gather, when the fft encounters (for ND > 1) even-numbered
% % % % % % dimensions, it dumps some complex conjugate pairs in a line at the edge
% % % % % % of the spectrum. For only odd dimensions, the each coeff sits happily with
% % % % % % its complex conjugate pair opposite it as a reflection through the very
% % % % % % centre of the spectrum, where the 0th freq component sits. If any
% % % % % % dimension is even, the fft dumps the extra coeffs in a line at the edge,
% % % % % % who sit opposite their complex conjugate pairs via a reflection with a
% % % % % % 0th freq component in the centre of that extra line (see below).
% % % % % %
% % % % % % The number of 2s should always equal the number of 0s.
% % % % % %
% % % % % % For a dimension with ODD number of elements, simply take the mask below
% % % % % % from inside the line, and for EVEN elements include the elements outside
% % % % % % of the lines. You essentially seem to trim the below to suit your needs
% % % % % % for each of these planes.
% % % % % %
% % % % % % MASK EXAMPLE:
% % % % % % EVEN n_elements case | ODD n_elements case
% % % % % %
% % % % % %         -ve X     +ve X
% % % % % %   2 | 0 0 0 0 2 2 2 2 2
% % % % % %   2 | 0 0 0 0 2 2 2 2 2 +ve Y
% % % % % %   2 | 0 0 0 0 2 2 2 2 2
% % % % % %   2 | 0 0 0 0 2 2 2 2 2
% % % % % %   1 | 0 0 0 0 1 2 2 2 2
% % % % % %   0 | 0 0 0 0 0 2 2 2 2
% % % % % %   0 | 0 0 0 0 0 2 2 2 2
% % % % % %   0 | 0 0 0 0 0 2 2 2 2 -ve Y
% % % % % %   0 | 0 0 0 0 0 2 2 2 2
% % % % % %   --+------------------
% % % % % %   1 | 0 0 0 0 1 2 2 2 2
% % % % % %
% % % % % %
% % % % % % Once you've set these 4 planes correctly, the rest is fairly
% % % % % % straightforward - simply set everything else above the middle horizontal
% % % % % % plane to be 2 (for +ve Z freqs) and set everything below (but not the
% % % % % % base) to be 0 (not -ve Zfreqs). You can of course fiddle this method to
% % % % % % consider different freq combinations, but tbh it's probably easier just
% % % % % % to permute you matrix to what you want and put it through this code :)
% % % % % %
% % % % %
% % % % % %%% Now the function itself!
% % % % % %
% % % % % % INPUTS: sz - size() of the desired mask, either 1, 2 or 3D.
% % % % % %
% % % % % % OUTPUTS: m - the mask itself, same size as input dimensions
% % % % % %
% % % % % % Note, we expect an fftshifted vector/matrix in each case.
% % % % % %
% % % % %
% % % % % function m = nph_hilbertmask(sz,varargin)
% % % % %
% % % % % % sz = size of the input matrix to be transformed.
% % % % % % negflag = 'neg' = some flag to denote that we need to create a negative
% % % % % % frequency-including mask on the 3rd dimension.
% % % % %
% % % % % % add scales, if supplied. This allows us to switch between +ve and -ve
% % % % % % masks for the 3rd dimension in the 3DST. Before, we only considered
% % % % % % positive z freqs.
% % % % % negflag = 0;
% % % % % if nargin == 2
% % % % %     if any(strcmpi(varargin{1},{'neg','negflag','-'}))
% % % % %         negflag = 1;
% % % % %     end
% % % % % end
% % % % %
% % % % %
% % % % %
% % % % % m = nan(sz); % ensure output is same size as specified
% % % % %
% % % % % sz(sz == 1) = []; % cope with 1x8, 8x1x5 etc.
% % % % %
% % % % % mid = fix(sz/2)+1; % find midpoint(s).
% % % % %
% % % % % switch length(sz)
% % % % %
% % % % %     %   1D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %     case 1 % 1D (straight from bob stockwell)
% % % % %
% % % % %         if isodd(sz)
% % % % %             m(1:mid-1)      = 0;
% % % % %             m(mid)          = 1;
% % % % %             m(mid+1:end)    = 2;
% % % % %         else
% % % % %             m(1)            = 1;
% % % % %             m(2:mid-1)      = 0;
% % % % %             m(mid)          = 1;
% % % % %             m(mid+1:end)    = 2;
% % % % %         end
% % % % %         %         m = [ones(1-rem(sz,2),1); 2*ones(fix((sz-1)/2),1); 1; zeros(fix((sz-1)/2),1)];
% % % % %
% % % % %         %   2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %     case 2 % 2D (NPH and NDS method)
% % % % %         % make the 2D mask featured in the preamble:
% % % % %         m = mask2D(sz);
% % % % %
% % % % %         %   3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %     case 3 % 3D (NPH method)
% % % % %
% % % % %         % First, select all +/-ve X,Y and +ve Z freqs:
% % % % %         % SWITCH CASE FOR NEGATIVE Z FREQS:
% % % % %         switch negflag
% % % % %             case 0 % +ve z freqs
% % % % %                 m(:,:,mid(3)+1:end) = 2;
% % % % %                 m(:,:,1:mid(3)-1)   = 0;
% % % % %             case 1 % -ve z freqs
% % % % %                 m(:,:,mid(3)+1:end) = 0;
% % % % %                 m(:,:,1:mid(3)-1)   = 2;
% % % % %         end
% % % % %
% % % % %         % Do the middle slice (always done regardless of odd/even
% % % % %         % dimensions)
% % % % %         m(:,:,mid(3)) = mask2D(sz([1 2]));
% % % % %
% % % % %         % Now, determine what extra sides you need depending on the
% % % % %         % odd/even dimensions:
% % % % %
% % % % %         % Base
% % % % %         if iseven(sz(3))
% % % % %             m(:,:,1) = mask2D(sz([1 2]));
% % % % %         end
% % % % %
% % % % %         % left hand side
% % % % %         if iseven(sz(2))
% % % % %             m(:,1,:) = mask2D(sz([1 3]));
% % % % %         end
% % % % %
% % % % %         % front
% % % % %         if iseven(sz(1))
% % % % %             m(1,:,:) = mask2D(sz([2 3]));
% % % % %         end
% % % % %
% % % % %
% % % % %
% % % % % end
% % % % %
% % % % % % Sanity check:
% % % % % if numel(m == 2) ~= numel(m == 0)
% % % % %     disp('Warning: Problem with the Hilbert mask.')
% % % % %     return
% % % % % end
% % % % %
% % % % %
% % % % % end
% % % % %
% % % % %
% % % % % % 2D MASK ==========================================================================
% % % % % function m2 = mask2D(dims)
% % % % %
% % % % % mid = fix(dims/2)+1;
% % % % %
% % % % % % the order matters!
% % % % % m2 = nan(dims);
% % % % %
% % % % % m2(:,mid(2):end)            = 2; % double half
% % % % % m2(1:mid(1),1:mid(2))       = 0; % zero bottom quarter
% % % % % m2(mid(1):end,1:mid(2)-1)   = 0; % zero top quarter
% % % % %
% % % % % m2(mid(1),mid(2))           = 1; % middle 0th freq
% % % % %
% % % % % % If both sides are odd, the 2D mask is finished now! yayyy!
% % % % %
% % % % % % Bottom and left hand side edges?
% % % % % switch double(iseven(dims(1))) / double(iseven(dims(2)))
% % % % %     case Inf % [1 0] odd, even
% % % % %         m2(1,:) = [zeros(1,mid(2)-1) 1 2*ones(1,mid(2)-1)];
% % % % %     case 0   % [0 1] even, odd
% % % % %         m2(:,1) = [zeros(1,mid(1)-1) 1 2*ones(1,mid(1)-1)]';
% % % % %     case 1   % [1 1] even, even
% % % % %         m2(1,:) = [1 zeros(1,mid(2)-2) 1 2*ones(1,mid(2)-2)];
% % % % %         m2(:,1) = [1 zeros(1,mid(1)-2) 1 2*ones(1,mid(1)-2)]';
% % % % % end
% % % % %
% % % % % % Sanity check:
% % % % % if numel(m2 == 2) ~= numel(m2 == 0)
% % % % %     disp('Warning: Problem with the Hilbert mask.')
% % % % %     return
% % % % % end
% % % % %
% % % % % % flipud(m2 - nph_hilb2(dims))
% % % % % % figure; imagesc(1:dims(1),1:dims(2),m2);
% % % % %
% % % % %
% % % % % end
% % % % % %==========================================================================


%==========================================================================
% ISEVEN and ISODD functions
function logikal = iseven(x)
logikal = mod(x,2)==0;
end
function logikal = isodd(x)
logikal = mod(x,2)==1;
end
%==========================================================================

% % % % % %
% % % % % % % CREATE DEFAULTS =========================================================
% % % % % % function Defaults = create_defaults(IN,type)
% % % % % %
% % % % % % switch type
% % % % % %     case 1 % 1DST
% % % % % %         % Default scales = 1:Nyquist for 1DST;
% % % % % %         Defaults.scales = 1:length(IN)-1;
% % % % % %         Defaults.point_spacing = 1;
% % % % % %         Defaults.c = 0.25;
% % % % % %
% % % % % %     case 2 % 2DST
% % % % % %         % Default scales = 1:N/3
% % % % % %         Defaults.scales = {...
% % % % % %             [-fix(size(IN,1)/3):1:-1 1:fix(size(IN,1)/3)], ...
% % % % % %             [-fix(size(IN,2)/3):1:-1 1:fix(size(IN,2)/3)]};
% % % % % %         Defaults.point_spacing = [1 1];
% % % % % %         Defaults.c = [0.25 0.25];
% % % % % %
% % % % % %     case 3 % 3DST
% % % % % %         % Default scales = 1:N/3
% % % % % %         Defaults.scales = {...
% % % % % %             [-fix(size(IN,1)/3):1:-1 1:fix(size(IN,1)/3)], ...
% % % % % %             [-fix(size(IN,2)/3):1:-1 1:fix(size(IN,2)/3)], ...
% % % % % %             [-fix(size(IN,3)/3):1:-1 1:fix(size(IN,3)/3)]};
% % % % % %         Defaults.point_spacing = [1 1 1];
% % % % % %         Defaults.c = [0.25 0.25 0.25];
% % % % % %
% % % % % % end
% % % % % %
% % % % % % end
% % % % % % %==========================================================================




% FFT WISDOM =============================================================
function fftwisdom = gatherfftwisdom(IN,ffttype)
% So we may be able to significantly improve the speed of the ifftn
% computation below by using the fft wisdom library:

% method = 'exhaustive';
method = 'patient';

switch class(IN)
    case 'double'
        wisdomtype = 'dwisdom';
    case 'single'
        wisdomtype = 'swisdom';
end

fftw(wisdomtype,''); % clear any existing wisdom
fftw('planner',method); % choose method

switch lower(ffttype)
    case {'fft' 'fftn'}
        test_fftn = fftn(IN); % test fftn
    case {'ifft' 'ifftn'}
        test_ifftn = ifftn(IN); % test ifftn
end

fftwisdom = fftw(wisdomtype); % get wisdom

fftw(wisdomtype,fftwisdom); % apply the wisdom


end


































