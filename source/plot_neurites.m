function [fig1,fig2] = plot_neurites( segments,varargin )
    % fig1 = plot_neurites( segments ) plots a montage of the individual segments
    % detected by "locate_neurites.m". The figure handle is stored in fig1
    %
    % fig1 = plot_neurites( segments,field ) sorts the segments by their
    % values specified by field. field is a string, and can be any one of the
    % fields in the "segments" structure that isn't a cell array. The
    % segments will be sorted by the values contained in this field in
    % increasing order. Default = no sorting. Set to [] to leave default
    %
    % fig1 = plot_neurites( ...,useBinarized ) uses the binarized (skeleton)
    % versions of the individual segments for the montage if set to true.
    % Default = false (uses the non-binarized versions)
    %
    % [fig1,fig2] = plot_neurites( ...,processed ) also plots the variosu processed
    % images in the "processed" structure, also computed from
    % "locate_neurites.m". fig2 is the figure handle
    %
    % written by Jordan Sorokin
    % 4/21/2018
    
    % check inputs
    sortSegs = false;
    if nargin > 1 && ~isempty( varargin{1} )
        field = varargin{1};
        validFields = fields( segments );
        cells = structfun( @iscell,segments );
        if ~any( ismember( validFields,field ) & ~cells )
            fprintf( '%s is not a valid name for sorting segments\n',varargin{1} );
            fprintf( 'resorting to no sorting\n' );
        else
            sortSegs = true;
            values = eval( ['segments.',field] );
        end
    end
    
    if nargin > 2 && ~isempty( varargin{2} )
        useBinarized = varargin{2};
    else
        useBinarized = false;
    end
    
    if nargin > 3 && ~isempty( varargin{3} )
        processed = varargin{3};
        plotProcessed = true;
    else
        processed = [];
        fig2 = [];
        plotProcessed = false;
    end
    
    warning( 'off' );
    
    %% plot the montage
    fig1 = figure('Name','Detected neurites','Position',[900 250 600 600]);
    
    if useBinarized
        images = segments.skeletons;
    else
        images = segments.images;
    end
    
    % sort if needed
    if sortSegs
        [~,idx] = sort( values,'ascend' );
        images = images(idx);
    end
    
    % make the segments equal size and plot
    bigImages = pad_neurites( images );
    montage( bigImages,'DisplayRange',[] ); axis image off; colormap(gray(256));
    
    if sortSegs
        title( ['Sorted by: ',field] );
    end
    %% plot the processed figs
    if plotProcessed
        fig2 = figure('Name','Processed original image','Position',[100 250 600 600],'Visible','off'); % processed fig
        ax(1) = subplot( 2,2,1 );
        imagesc( processed.pooledImg,[0,max(processed.pooledImg(:))] ); axis image off
        title( sprintf( 'gabor filter (%s pooling)',processed.params.pooling ) );

        ax(2) = subplot( 2,2,2 );
        imagesc( processed.cleanedImg,[0,max(processed.cleanedImg(:))] ); axis image off
        title( 'cleaned image' );

        ax(3) = subplot( 2,2,3 );
        imagesc( processed.skeletonImg,[0,1] ); axis image off
        title( 'detected ridges' );

        ax(4) = subplot( 2,2,4 );
        imagesc( processed.labeledImg,[0,segments.nSegs] ); axis image off
        title( 'detected components' );
        
        linkaxes( ax,'xy' );
        colormap( 'gray' );
    
        fig2.Visible = 'on';
    end
    warning( 'on' );
end
    
            