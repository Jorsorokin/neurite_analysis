function [segments,processed] = locate_neurites( img,degRange,varargin )
    % [segments,processed] = locate_neurites( img,degRange,varargin )
    %
    % locates filament-like objects in a 2D grayscale image "img"
    %
    % ASSUMPTIONS: neurites will be roughly oriented along some major axis,
    % (thus are not circular objects). The gabor filters do a great job at
    % filtering out blob-like structures, while preserving filaments
    %
    % Inputs:
    %   img:
    %       2D image, a projection of a possibly 3D stack
    %   
    %   degRange:
    %       range of angles for the gabor filter bank (in degrees...90* =
    %       horizontal, 0* = vertical)
    %
    %   (wavelength):
    %       the wavelength for the gabor. smaller = courser gabor filter
    %       (smoother but potentially blotchier filtering). (Default=10)
    %
    %   (gaborFrequency):
    %       number > 0 indicating the sptaial frequency of the gabor.
    %       Larger ( > 1 ) results in gabors that pick out higher
    %       frequency components in the image. (Default=1.5)
    %
    %   (aspectRatio):
    %       defines the major-to-minor axis ratio of the gabor filter. An
    %       aspect ratio of 1 means both the major and minor axes of the
    %       filter are equal. This will pick out edges equally in the x-
    %       and y- directions. For objects with a long spatial scale in one
    %       direction and less so in the other (like AIS), use an aspect
    %       ratio < 1, which elongates the gabor along the major axis (the
    %       axis defined by the angle of the filter). (Default = 0.5)
    %   
    %   (minPixels):
    %       minimum # of pixels needed for an object to be kept
    %       (default=10)
    %
    %   (alpha):
    %       level for changing the thresholding (smaller = keep dimmer
    %       objects). Default = 1
    %
    %   (gaussLength):
    %       the radius of the gaussian filter for ridge detection
    %       (larger=smoother but reduces spatial resolution). (default=15)
    %
    %   (gaussSigma):
    %       the SD of the gaussian filter (larger=lower frequency filter response).
    %       (default=3)
    %   
    %   (pooling):
    %       the method used to cumulate across gabor filters. Can be one of
    %       the following:
    %           'avg' - takes the average pixel intensity across all
    %                   filtered images
    %           'wAvg' - takes a weight exponential average by finding the 
    %                    filter that maximizes the pixel intensity for
    %                    each pixel individually, then weights the other
    %                    filters by a gaussian of size "tau" (see below)
    %           'max' - takes the maximum value for each filter / pixel
    %           'med' - takes the median value for each pixel / filter
    %       (default='max')
    %
    %   (tau) - sets the size of the gaussian for the weighted average
    %           across filters. For instance, if you have 10 filtered
    %           versions of "img", and for pixel j you find the 5th
    %           filtered image produced the brightest intensity, and tau is
    %           = 0.3, then the pixel j will be averaged across the 4th, 5th,
    %           and 6th filtered version as: 
    %
    %                   pixel_j = 0.1664*pixel_j4 + 0.6672*pixel_j5 + 0.1664*pixel_j6
    %           
    %           where pixel_jN means the value of pixel_j on the Nth filtered image 
    %
    %           To understand the coefficients used here, "tau" is percent of the total 
    %           # of gabor filters that we have. So here, we have 10
    %           gabor filters, and tau = 0.3, so our averaging is over
    %           10*0.3 = 3 filters. Then, the coefficients are created via:
    %               
    %                   "gausswin( tau )" and dividided by it's sum. 
    %   
    %           (Default=0.25)
    %
    % Outputs:
    %   segments:
    %       a structure with the following fields
    %           nSegs           -   number describing the # of AIS detected
    %
    %           pixels          -   a cell array, with each element being an N x 2 vector describing 
    %                               the X/Y pixel locations each object
    %
    %           nPixels         -   N x 2 vector describing the total # of X and Y
    %                               pixels for each image (row)
    %
    %           area            -   vector describing the # of pixels total for each segment
    %
    %           gaborAngle      -   cell array with each cell containing
    %                               the gabor angle that produced the
    %                               maximum response for each pixel of that segment
    %                               NOTE: 90* is subtracted
    %
    %           meanAngle       -   vector describing the mean angle in degrees
    %                               of each segment (i.e. gabor angle with maximum
    %                               response). NOTE: 90* is subtracted
    %
    %           curvature       -   cell array describing the curvature at each pixel for each 
    %                               individual segment
    %
    %           meanCurvature   -   vector describing the mean amount of curvature
    %
    %           eccentricity    -   vector describing the eccentricity (proxy for curvature)
    %
    %           orientation     -   vector with the estimated orientation of
    %                               the segments (redundant with "meanAngle")
    %
    %           majorAxisLength -   # of pixels of the major axis 
    %
    %           minorAxisLength -   # of pixels of the minor axis
    %
    %           axisRatio       -   (majorAxisLength - minorAxisLength) / nPixels
    %
    %           linearSlope     -   the slope of the best-fit line 
    %
    %           residuals       -   cell array of residuals to the fit
    %
    %           mse             -   the mean squared error of the residuals
    %
    %           centerOfMass    -   center of mass of the fluorescence of
    %                               the pixels in the cleaned segment
    %                               caluclated as: 
    %                                       sum( img,axis ) / sum( img(:) ) * [1:length( axis )]'
    %                               for each axis x and y
    %
    %           cmRatio         -   the relative location of the center of
    %                               mass on each segment calculated as:
    %                                   sum( centerOfMass / nPixels ) / 2
    %                               cmRatios at the extremes (0, 1)
    %                               indicate most fluorescence is at one
    %                               end of the segment
    %
    %           images          -   cell array of small snippet of the cleaned image created
    %                               by taking a bounding box around each segment
    %
    %           skeletons       -   same as above, but just the binary "ridges" of
    %                               each segment, which is used for the above
    %                               quantifications 
    %
    %   processed:
    %       another structure containing the results of the various
    %       processing steps on the input image:
    %           poooledImg - the result from pooling the gabor-filtered
    %                        images for the different gabors based on the
    %                        "pooling" method
    %           cleanedImg - a cleaned version of the input image (used for
    %                        detecting the skeletons of the segments)
    %           skeletonImg - a binarized version of cleanedImg, with
    %                         "ridges" for the detected segments
    %           labeledImg - skeletonImg, color coded by each detected
    %                        segment via "bwconncomp()"
    %           params - a structure with the values for the optional inputs
    %
    % Written by Jordan Sorokin
    % 4/21/2018 - finished 1st draft of script
    % 5/2/2018 - modified ridge detection and cleaned up segments
    % 5/2/2018 - added "centerOfMass" and "cmRatio" 
    
    %%
    
    warning( 'off' );
    
    % parse our inputs
    options = parse_inputs( varargin );
    
    % get size of image
    [n,m] = size( img );
    
    % create our filter banks and filter the image
    fprintf( 'computing responses to gabor filters...\n' );
    nGabors = numel( degRange );
    gaborFilters = gabor( options.wavelength,degRange,...
                          'SpatialFrequencyBandwidth',options.gaborFrequency,...
                          'SpatialAspectRatio',options.aspectRatio );
                      
    gaborFiltered = imgaborfilt( img,gaborFilters );
    [pooledImg,idx] = max( gaborFiltered,[],3 );
    
    % pool the filtered images using the specified method
    switch options.pooling
        case 'avg'
            pooledImg = mean( gaborFiltered,3 ); 
        case 'wAvg'
            tau = ceil( nGabors * options.tau );
            if mod( tau,2 ) ~= 0 
                tau = tau + 1; % change to even #
            end
            coeffs = gausswin( tau ); 
            coeffs = coeffs / sum( coeffs ); 
            nCoeffs = tau / 2; % for extracting pixels
            
            % for each pixel, find which filtered version of the image gave
            % the brightest pixel, then take a weighted average across
            % neighboring filtered versions. This is slow, but will better
            % preserve local structure vs. the "avg" method
            N = n*m;
            for i = 1:N
                bestFilter = idx(i);
                prevInds = max( bestFilter-nCoeffs,1 );
                postInds = min( bestFilter+nCoeffs,tau );
                [j,k] = ind2sub( [n,m],i );
                pooledImg(i) = mean( gaborFiltered(j,k,prevInds:postInds),3 );
            end
        case 'med'
            pooledImg = median( gaborFiltered,3 );
        case 'max'
    end
    
    % convert to binary image
    pooledImgBW = threshold( pooledImg,options.alpha );
    
    % now use the pooled thresholded image to denoise the original 
    % image. This produces a final thresholded image where we've eliminated 
    % segments that do not match our gabor filters
    cleanImg = img .* pooledImgBW; % cleans up the original to remove noise and non-aligned segments
    cleanImg = imclose( cleanImg,strel('square',3) ); % helps close line segments together
    cleanImgBW = cleanImg > quantile( cleanImg(:),0.9 ); % thresholds at 90% quantile (softer than 1 SD)
    
    % now compute the ridges in the clean image, using local Hessian
    % matrices to compute ridges in the intensity profile of the objects
    fprintf( 'skeletonizing the image...\n' );
    skeletonImg = ridgefilt( cleanImg,options.gaussLength,options.gaussSigma,0.85 ); % 85 % quantile works well
    skeletonImg = skeletonImg & cleanImgBW;
    [~,skeletonImg] = edgelink( skeletonImg ); % cleans up end points to remove blobs
    skeletonImg = bwareaopen( skeletonImg,options.minPixels ); % removes small objects
    
    % detect connected components
    fprintf( 'detecting and quantifying segments...\n' );
    components = bwconncomp( skeletonImg,8 ); % 8 connectivity for diagonal detection
    properties = regionprops( components,{'Area','MajorAxisLength','MinorAxisLength',...
                                          'SubarrayIdx','Eccentricity','Orientation','Image'});
    
    % clean up the connected components by removing circular (blobs)
    % without an orientation to them
    artifacts = [properties.Eccentricity] < 0.75; % low eccentricity = circular structures
    artifacts = artifacts | [properties.MajorAxisLength] ./ [properties.MinorAxisLength] <= 1.5;
    components.PixelIdxList(artifacts) = [];
    components.NumObjects = nnz( ~artifacts );
    properties(artifacts) = [];
    
    % re-binarize the image to get the final skeletonized version
    labeledImg = labelmatrix( components );
    skeletonImg = labeledImg ~= 0;
    
    % extract the segments for the raw image and binarized images, as well
    % as the various properties that we've extracted
    nSegs = components.NumObjects;
    pixelInds = components.PixelIdxList;
    [rows,columns] = cellfun( @(x)(ind2sub( [n,m],x )),pixelInds,'un',0 );
    pixels = cellfun( @(x,y)([x,y]),rows,columns,'un',0 );
    
    angles = cellfun( @(x)(degRange(idx(x))),pixelInds,'un',0 ); % gabor angle response
    curvature = cellfun( @(x)(LineCurvature2D( x )),pixels,'un',0 ); % curvature 
    [fitParams,residuals,mse] = cellfun( @(x)(linear_residuals( x )),pixels,'un',0 ); % slope and residuals of linear fit

    % pull out the segment from the clean image
    subfig = cell( 1,nSegs );
    cm = zeros( nSegs,2 );
    nPixels = zeros( nSegs,2 );
    for c = 1:nSegs
        i = pixels{c}(:,1); j = pixels{c}(:,2);
        xInds = min(i):max(i);
        yInds = min(j):max(j);
        subfig{c} = cleanImg( xInds,yInds ) .* imdilate( properties(c).Image,strel('disk',4) ); % removes pixels not part of segment
        tempImg = subfig{c};
        nPixels(c,:) = [numel( yInds ), numel( xInds )];
        
        % center of mass 
        xB = sum( tempImg,1 );
        yB = sum( tempImg,2 );
        cm(c,:) = [xB/sum(xB) * (1:nPixels(c,1))', (1:nPixels(c,2)) * yB/sum(yB)];
    end
       
    % store the results
    segments = struct;
    segments.nSegs = nSegs;
    segments.pixels = pixels;
    segments.nPixels = nPixels;
    segments.area = [properties.Area];
    segments.gaborAngle = cellfun( @(x)(x - 90),angles,'un',0 );
    segments.meanAngle =  cellfun( @nanmean,segments.gaborAngle );
    segments.curvature = curvature;
    segments.meanCurvature = cellfun( @nanmean,curvature );
    segments.orientation = [properties.Orientation];
    segments.eccentricity = [properties.Eccentricity];
    segments.majorAxisLength = [properties.MajorAxisLength];
    segments.minorAxisLength = [properties.MinorAxisLength];
    segments.axisRatio = (segments.majorAxisLength - segments.minorAxisLength ) ./ segments.area;
    segments.linearSlope = cell2mat( fitParams );
    segments.residuals = residuals;
    segments.mse = cell2mat( mse );
    segments.centerOfMass = cm;
    segments.cmRatio = mean( segments.centerOfMass ./ segments.nPixels,2 );
    segments.images = subfig;
    segments.skeletons = {properties.Image};
    
    processed = struct;
    processed.pooledImg = gather( pooledImg );
    processed.cleanedImg = cleanImg;
    processed.skeletonImg = skeletonImg;
    processed.labeledImg = labeledImg;
    processed.params = options;
    processed.params.gaborAngles = degRange;

    warning( 'on' );
    %% FUNCTIONS 
    function inputs = parse_inputs( optionalInputs )
        names = {'wavelength','gaborFrequency','aspectRatio','minPixels',...
            'pooling','tau','alpha','gaussLength','gaussSigma'};
        defaults = {10,1.5,0.5,10,'max',0.5,1,15,3};
        
        p = inputParser();
        for q = 1:numel( names )
            p.addParameter( names{q},defaults{q} );
        end
        
        p.parse( optionalInputs{:} );
        inputs = p.Results;
    end

    function imageBW = threshold( img,alpha )
        % thresholds the image to produce a black-white binary image
        threshLev = mean( img(:) ) + std( img(:) );
        imageBW = img > threshLev*alpha;
    end

    function [P,residuals,mse] = linear_residuals( pixels )
        % fit 1D polynomial to the pixels (a line)
        x = flipud( pixels(:,1) ); y = pixels(:,2); % flipped x to match original image
        M = numel( x );
        P = polyfit( x,y,1 );
        
        % get the fitted line, then compute residuals along the pixels
        curve = polyval( P,x );
        residuals = y - curve;
        
        % now standardize the residuals
        X = [ones( M,1 ),x];
        H = X*((X'*X)\X');
        h = diag( H );
        sigma = 1/(M-2) * sum( residuals.^2 );
        residuals = residuals ./ sqrt( sigma * (1-h) );
        
        % compute the mean squared residual
        mse = mean( residuals.^2 );  
        
        % extract just the slope
        P = P(1);
    end
end
