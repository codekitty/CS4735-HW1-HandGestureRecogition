function identify_password(filename)

% clean workspace
clc;
% close all;
verbos = true;

% init variables
pos = [];

% If left unspecified, ask user to supply a movie clip in '.mov' format
if ~exist('filename', 'var')
    [filename, pathname, ~] = uigetfile('*.mov', 'Load Passowrd Video');
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
    
    % Detect Skin
%     Cg = vidFrame(:,:,1);
%     Cb = vidFrame(:,:,2);
%     Cr = vidFrame(:,:,3);
%     binFrame = Cg>=120 & Cb>=130 & Cb<=255 & Cr>=130 & Cr<=255;
 
    % Segment image
    [thresh, em] = graythresh(vidFrame);
        
    [binFrame, contours, regions] = segment_image(vidFrame, thresh);
    
    % locate hand in the many continous regions
    if (length(contours)>1)
    
        % check if there is more than one large region - go tougher on thresh
        while length(contours{1})/2 < length(contours{2})
            thresh = thresh*1.05; % increase by five percent
            [binFrame, contours, regions] = segment_image(vidFrame, thresh);
            if (thresh>.9)
                break;
            end
            if verbos
                fprintf('\nwarning: increased thresh to %g', thresh);
            end
        end

        % check for hand in the two largest continous regions
        handidx =1;
        if ~isHand(contours{1}, binFrame, pos) && isHand(contours{2}, binFrame, pos) 
            if verbos
                fprintf('\nwarning: selected second largest region');
            end
            handidx = 2;
        end

        contHand = contours{handidx};
%     regHand  = regions{handidx}

    else
        contHand = contours{1};
        regHand  = regions;
    end
    
    % Where: center of mass
    pos(end+1, :) = round(mean(contHand));
    
    % what: use convex hull
    ch = convhull(contHand);
    
    % search for defects
    chpoints = conthand(ch, :);
    defect_score = [];
    for p=2:length(ch)
        d = norm(contHand(ch(p-1), :)-contHand(ch(p), :));
        [k, v] = convhulln([contHand(ch(p-1), :); contHand(ch(p), :); ]);
        defect_score(p-1) = v/d;
    end
    didx = sort(defect_score);
    
    if (verbos)
        [i,j,s] = find(binFrame);
        scatter(j, i, 100, '.b');    
        scatter(pos(:,2), pos(:,1), 400, '.r');
        scatter(contHand(:,2),contHand(:,1), 100, '.k');
        plot(contHand(ch,2), contHand(ch,1), '-r');
        title(['em ' num2str(em)]);
    end
    
     pause(1/(vidObj.FrameRate*100));
     hold off;
end

end

% check if contour is the same object we've been tracking lately
function hand=isHand(contour, img, past_positions)
% current logic: check if the new center of mass is within ten percent of
% the median of the last ten positions.
if size(past_positions,1)<11
    hand=true;
    return;
end

hand = norm(mean(contour)-median(past_positions(end-10:end, :))) < length(img)/10;
end

function [binFrame, contours, regions]=segment_image(img, thresh)
    % convert to pos\neg image
    binFrame = im2bw(img, thresh);

    % continous regions - contours
    [contours,regions] = bwboundaries(binFrame,'noholes');
    
    % sort by size
    [sorted, idx] = sort(cellfun(@(x)length(x),contours), 'descend'); 
    contours = contours(idx);
%     regions = regions(idx);
end