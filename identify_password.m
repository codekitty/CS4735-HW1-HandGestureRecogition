function identify_password(filename)

% params
verbos = true;
pink_skin = true;

% clean workspace
clc;
close all;

% init variables
pos = [];

% If left unspecified, ask user to supply a movie clip in '.mov' format
if ~exist('filename', 'var')
    [filename, pathname, ~] = uigetfile('*.mov', 'Load Password Video');
    if isequal(filename,0) || isequal(pathname,0)
        return;
    end
end

% open video
vidObj = VideoReader(filename);

% while video still has frames
while hasFrame(vidObj)
    % read new frame
    vidFrame = readFrame(vidObj);    
    imshow(vidFrame);
    hold on;
     
    % Segment image
    [thresh, em] = graythresh(vidFrame);
     
    % convert to pos\neg image
    binFrame = rgb2bin(vidFrame, thresh, pink_skin);
    [contours, regions] = segment_image(binFrame);
    
    % locate hand in the many continous regions
    if (length(contours)>1)
    
        % check if there is more than one large region - go tougher on thresh
        while length(contours{1})*.8 < length(contours{2})
            thresh = thresh*1.05; % increase by five percent
            binFrame = rgb2bin(vidFrame, thresh, pink_skin);
            [contours, regions] = segment_image(binFrame);
            if (thresh>.9)
                break;
            end
            if verbos
                fprintf('\nwarning: increased thresh to %g', thresh);
            end
        end

        % check for hand in the two largest continous regions
        handidx = 1;
        reg1 = getInds(regions==1);
        reg2 = getInds(regions==2);
        if length(reg1)/2 < length(reg2) && ...
            ~isCloseToRecent(reg1, binFrame, pos) && isCloseToRecent(reg2, binFrame, pos) 
            if verbos
                fprintf('\nwarning: selected second largest region');
            end
            handidx = 2;
        end

        contHand = contours{handidx};
        binHand  = (regions==handidx);
    else
        contHand = contours{1};
        binHand  = (regions==1);
    end
    
    % Where: center of mass
    [i,j,~] = find(binHand);
    pos(end+1, :) = round(mean([i, j]));
    
    % what: use convex hull
    ch = convhull(contHand);
    
    % search for convexity defects
    defects = [];
    for p=3:length(ch)
        if abs(ch(p-1) - ch(p)) < 35;
            continue; % must have at some pixels inbetween to be considered.
        end
        d = norm(contHand(ch(p-1), :)-contHand(ch(p), :));
        from=min(ch(p-1), ch(p));
        to=max(ch(p-1), ch(p));
        nuk = contHand(from:to, :);
        nukch = convhull(nuk);
        [nn, dists] = knnsearch(nuk(1, :), nuk(nukch(1:(end-2)), :));
        if any(dists >= d)  
            [~,maxdistsid] = max(dists);
            defects(end+1,:) = nuk(nukch(maxdistsid), :);
            
        end
    end
    
    % estimate a circle around palm
    if ~isempty(defects)
        [~, r] = knnsearch(pos(end, :),defects); 
        palm_r = max(r);
    end
        
    % plot the frame with identified objects and metrics
    if (verbos)
%         frameInds = getInds(binFrame);
%         scatter(frameInds(:,2), frameInds(:,1), 100, '.b');    
        scatter(pos(:,2), pos(:,1), 400, '.r');
        scatter(contHand(:,2),contHand(:,1), 100, '.w');
        if ~isempty(defects)
            scatter(defects(:, 2), defects(:, 1), 400,'.g');
            plot_circle(pos(end, :),palm_r, '-y');
        end
        plot(contHand(ch,2), contHand(ch,1), '-r');
        title(['em ' num2str(em)]);
    end
    
    pause(1/(vidObj.FrameRate*100));
    hold off;
end

% == analyse for gesture pass key ==

% examine the vectors of directions
dirs = pos(2:end, :)-pos(1:(end-1), :);

% median filter for reducing noise
spos = smooth(pos, 20);
scatter(spos(:,2), spos(:,1), 400, '.m');

end

% check if contour is the same object we've been tracking lately
function hand=isCloseToRecent(contour, img, past_positions)
% current logic: check if the new center of mass is within ten percent of
% the median of the last ten positions.
if size(past_positions,1)<11
    hand=true;
    return;
end

hand = norm(mean(contour)-median(past_positions(end-10:end, :))) < length(img)/20;
end

% segments an image to continous regions. regions\contours are returned
% sorted in descending length of perimeter.
function [contours, regions]=segment_image(binFrame)

    % continous regions - contours
    [contours,orig_regions] = bwboundaries(binFrame,'noholes');
    
    % sort by size
    [sorted, idx] = sort(cellfun(@(x)length(x),contours), 'descend'); 
    contours = contours(idx);
    
    % sort region numbers
    regions= orig_regions;
    for i=1:numel(idx)
        regions(orig_regions==idx(i)) = i;
    end
end

% converts an rgb image to a binary matrix where 1 identifies pixels of a
% hand and zeros are background
function binFrame = rgb2bin(img, thresh, pink_skin)
    binFrame = im2bw(img, thresh);
    
    if (pink_skin)
        % Detect Skin
        Cr = img(:,:,1);
        Cg = img(:,:,2);
        Cb = img(:,:,3);
        binFrame = binFrame & Cr>(Cg+10) & Cr>(Cb+10) & Cb<=230;
    end
end

% plots a circle with a given center and radius on top of image
function plot_circle(center, radius, style)
    n = 100;
    t=linspace(0,2*pi,n);
    r=ones(1,n)*radius;
    [x,y] = pol2cart(t,r);
    x=x+center(1);
    y=y+center(2);
    plot(y,x,style);
end

% returns a two-column matrix with indices corresponding to the indices
% where the given matrix is positive
function inds=getInds(mat)
    [i,j,~] = find(mat);
    inds = [i,j];
end

