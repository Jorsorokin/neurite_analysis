function bigImages = pad_neurites( images )
    % bigImages = pad_neurites( images ) pads the individual neurites
    % contained in the cell-array "images" with zeros so that all images
    % are of equal size. The actual neurites will be centered within each
    % padded image. useful for montage display / PCA
    %
    % Written by Jordan Sorokin, 5/1/18
    
    type = class( images{1} );
    K = numel( images );
    [n,m] = cellfun( @size,images);
    nMax = max( n );
    mMax = max( m );
    rowPad = nMax - n;
    colPad = mMax - m; 

    top = floor( rowPad/2 );
    bottom = top + (mod( rowPad,2 ) ~= 0);
    left = floor( colPad/2 );
    right = left + (mod( colPad,2 ) ~= 0);
    bigImages = zeros( nMax,mMax,1,K,type );

    for k = 1:K
        img = images{k};
        leftPad = zeros( n(k),left(k),type );
        rightPad = zeros( n(k),right(k),type );
        topPad = zeros( top(k),m(k)+colPad(k),type );
        bottomPad = zeros( bottom(k),m(k)+colPad(k),type );

        newImg =       [ topPad;
                [leftPad, img, rightPad];
                       bottomPad ];

        bigImages(:,:,1,k) = newImg;
    end
end