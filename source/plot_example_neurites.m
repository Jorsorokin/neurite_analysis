function fig = plot_example_neurites( segments,sortFields,quantiles,varargin )
    % fig = plot_example_neurites( segments,sortFields,quantiles )
    % plots example neurites sorted by potentially multiple fields. Each
    % field is 1 row, each column is the representative quantile
    %
    % fig = plot_sexample_neurites( ...,useBinarized ) plots the binarized
    % versions of the segments instead of the cleaned non-binarized ones.
    % useBinarized is a boolean flag (default = false)
    %
    % % EX
    % % plot an example neurite for the 10, 50, and 90 % quantiles 
    % % sorted by "eccentricity" and "axisRatio". 
    % plot_neurites( segments,{'eccentricity','axisRatio'},[0.1, 0.5, 0.9] )
    % 
    % By Jordan Sorokin, 4/25/18
    
    % check if the fields provided are in the structure "segments"
    allFields = fields( segments );
    cells = structfun( @iscell,segments );
    allFields = allFields( ~cells );
    validFields = contains( allFields,sortFields );
    
    if ~any( validFields )
        disp( 'No valid fields for sorting provided' );
        fig = [];
        return
    end
    
    % set up the figure
    nRows = nnz( validFields );
    nCols = numel( quantiles );
    height = 0.9 / (nRows + 0.08*nRows);
    width = 0.9 / (nCols + 0.08*nCols);
    yPos = linspace( 0.95-height,0.05,nRows );
    xPos = linspace( 0.05,0.95-width,nCols );
    fig = figure( 'Visible','off','Position',[300 200 800*(min( nCols/nRows,1 )) 800*(min( nRows/nCols,1 ))] );
    
    if nargin > 3 && ~isempty( varargin{1} )
        useBinarized = varargin{1};
    else 
        useBinarized = false;
    end
    
    % loop over valid fields, plot each as a row
    allFields = allFields( validFields );
    if useBinarized
        images = segments.skeletons;
    else
        images = segments.images;
    end
    
    % pull out images so we can pad later
    examples = cell(1,nRows*nCols);
    counter = 1;
    for i = 1:nRows
        evalc( ['vals = segments.',allFields{i}] );
        for q = 1:nCols
            [~,idx] = min( (vals - quantile( vals,quantiles(q) )).^2 );
            examples{counter} = images{idx};
            counter = counter + 1;
        end
    end
    
    % pad the images and plot
    bigImages = squeeze( pad_neurites( examples ) );
    quantileCounter = 1;
    fieldCounter = 1;
    for i = 1:nRows*nCols
        ax = subplot( 'Position',[xPos(quantileCounter),yPos(fieldCounter),width,height] );
        imagesc( bigImages(:,:,i) );
        title( sprintf( '%s %0.2f quantile',allFields{fieldCounter},quantiles(quantileCounter) ) );
        set( ax,'fontsize',8); axis image off square
        if mod( i,nCols ) == 0
            quantileCounter = 1;
            fieldCounter = fieldCounter + 1;
        else
            quantileCounter = quantileCounter + 1;
        end
    end
    
    fig.Visible = 'on';
    colormap( 'bone' );
end