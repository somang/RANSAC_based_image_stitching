%clc;clear;
%run('C:\Users\imdc\Documents\MATLAB\SIFT\vlfeat-0.9.18\toolbox\vl_setup.m')
%run('/home/somang/Documents/SIFT/vlfeat-0.9.18/toolbox/vl_setup')

%% PART 1.
orig_left = imread('parliament-left.jpg');
orig_right = imread('parliament-right.jpg');
left = im2single(rgb2gray(orig_left));
right = im2single(rgb2gray(orig_right));

%% PART 2.
% Given:
% pix1 = [x1; y1];
% pix2 = [x2; y2];
% T = [t1 t2; t3 t4]; % 2 by 2
% c = [c1; c2] % 2 by 1
% Such that,
% pix2 = T*pix1 + c
% 
% Then, we can compute each element as
% [x2; y2] = [t1 t2; t3 t4]*[x1; y1] + [c1; c2]
%          = [x1t1 + y1t2 + c1; x1t3 + y1t3 + c2]
% 
% Where we have, 
% x2 = x1t1 + y1t2 + c1 AND y2 = x1t3 + y1t3 + c2
% 
% To compute the same result in Ax = b format, where x is a vector consists of unknowns.
% We can have;

% x = [t1; t2; t3; t4; c1; c2];
% A = [x1 y1 0 0 1 0; 0 0 x1 y1 0 1];
% b = [x2; y2];
% 
% Where once we compute Ax = b gives the same results.
% If we have multiple points, then we add more rows in our matrix A by
% adding other pixels coordinates rows below.
% A = [x1 y1 0 0 1 0; 0 0 x1 y1 0 1; (...) xn yn 0 0 1 0; 0 0 xn yn 0 1];

%% PART 3.
% Three pairs at least, since it needs to compare for the plane.

%% PART 4.
A = [197 2014 0 0 1 0; ...
     0 0 197 2014 0 1; ...
     1271 2183 0 0 1 0; ...
     0 0 1271 2183 0 1; ...
     1580 2010 0 0 1 0; ...
     0 0 1580 2010 0 1];
b = [97; 573; 1163; 685; 1455; 500];
x = A\b; % Solve for x

%% PART 5.
% There might be a way using xdata and ydata or using c1 c2 calculation...
intersection = 1160; % manual padding for better looking.
T = maketform('affine', [x(1) x(2) 0; x(3) x(4) 0; x(5) x(6) 1]);

%% PART 6.
left_t = imtransform(left, T);
composition = [ left_t zeros(size(left_t,1),intersection) ]; % put left ONTO right.
right_t = [ zeros(200, size(right,2)); right ]; 
right_window = size(right,2) - intersection; % intersecting window size.
left_window = size(composition,2) - intersection; % again.
intersect_pt = left_window-right_window-1; % Where it intersects.
composition(1:2600, left_window:end) = ...
    right_t(1:2600,right_window:end); % ONTO
composition(1:2600,intersect_pt:left_window) = ...
    composition(1:2600,intersect_pt:left_window)/2 + ...
    right_t(1:2600,1:right_window+2)/2;
 
%% PART 7.
r1 = orig_left(:,:,1); g1 = orig_left(:,:,2); b1 = orig_left(:,:,3);
r2 = orig_right(:,:,1); g2 = orig_right(:,:,2); b2 = orig_right(:,:,3);

left_t = imtransform(r1, T);
img_comp_r = [left_t zeros(size(left_t,1),intersection)];
right_t = [zeros(200, size(r2,2)); r2];
img_comp_r(1:2600,left_window:end) = ...
    right_t(1:2600,right_window:end);
img_comp_r(1:2600,intersect_pt:left_window) = ...
    img_comp_r(1:2600,intersect_pt:left_window)/2 + ...
    right_t(1:2600,1:right_window+2)/2;

left_t = imtransform(g1, T);
img_comp_g = [left_t zeros(size(left_t,1),intersection)];
right_t = [zeros(200, size(g2,2)); g2];
img_comp_g(1:2600,left_window:end) = ...
    right_t(1:2600,right_window:end);
img_comp_g(1:2600,intersect_pt:left_window) = ...
    img_comp_g(1:2600,intersect_pt:left_window)/2 + ...
    right_t(1:2600,1:right_window+2)/2;

left_t = imtransform(b1, T);
img_comp_b = [left_t zeros(size(left_t,1),intersection)];
right_t = [zeros(200, size(b2,2)); b2];
img_comp_b(1:2600,left_window:end) = ...
    right_t(1:2600,right_window:end);
img_comp_b(1:2600,intersect_pt:left_window) = ...
    img_comp_b(1:2600,intersect_pt:left_window)/2 + ...
    right_t(1:2600,1:right_window+2)/2;

imcolor = cat(3,img_comp_r,img_comp_g,img_comp_b);
figure, imshow(imcolor), pause;


%% PART 8,9.
% f(1:2) -> centre x,y
[f1,d1] = vl_sift(left);
[f2,d2] = vl_sift(right);

%% PART 10.
% matches f1, f2 and scores.
[matches, scores] = vl_ubcmatch(d1, d2);

%% PART 11.
filtered_match = zeros(size(matches,2), 2);
index = 1;
for i = 1:size(matches,2)
    if(scores(i) < 100) % Euclidean distances between SIFT descriptors are less than 100.
        filtered_match(index,1) = matches(1,i);
        filtered_match(index,2) = matches(2,i);
        index = index + 1;
    end
end
filtered_match = filtered_match(1:index-1,:);

%% PART 12.
% Implement RANSAC to estimate an ane transformation mapping one image onto the other. 
% Use the minimum number of pairwise matches to estimate the ane transformation. 
% Since you are using the minimum number of pairwise points, 
% the transformation can be estimated using an inverse transformation rather than least-squares.

best_fit_est = [0 0 0 0]; % will contain, the selected ransac pts with the most number of inliners contained.
selected_inliners = [];

for N = 1:100 % Let's repeat 100 times.
    ransac_pts = randperm(size(filtered_match,1),3); %pick three unique random numbers.
    % Get the coordinates.
    p1f1 = f1(1:2, filtered_match(ransac_pts(1),1));
    p1f2 = f2(1:2, filtered_match(ransac_pts(1),2));
    p2f1 = f1(1:2, filtered_match(ransac_pts(2),1));
    p2f2 = f2(1:2, filtered_match(ransac_pts(2),2));
    p3f1 = f1(1:2, filtered_match(ransac_pts(3),1));
    p3f2 = f2(1:2, filtered_match(ransac_pts(3),2));
    % Estimate the affine transformation
    A = [p1f1(2),p1f1(1),0,0,1,0; ...
        0,0,p1f1(2),p1f1(1),0,1; ...
        p2f1(2),p2f1(1),0,0,1,0; ...
        0,0,p2f1(2),p2f1(1),0,1; ...
        p3f1(2),p3f1(1),0,0,1,0; ...
        0,0,p3f1(2),p3f1(1),0,1];
    b = [p1f2(2); p1f2(1); p2f2(2); p2f2(1); p3f2(2); p3f2(1)];
    x = A\b;
    % Now, we have our unknowns.
    T = [x(1) x(2); x(3) x(4)];
    c = [x(5); x(6)];
    
    % check all other pts are inliners : i.e. has less distances than
    % threshold p.
    p = 0.025;
    inliners = []; % Contains which index of filtered_match has inliners.
    for i = 1 : size(filtered_match,1)
        if (ismember(i,ransac_pts) == 0) % if the point is not picked before
            tmp_pix = f1(1:2, filtered_match(i,1));
            pix = [tmp_pix(2); tmp_pix(1)];
            pix_map = T*pix+c; % Where does this map to?
            
            tmp_orig = f2(1:2, filtered_match(i,2)); % what was the real match?
            
            % Calculate Euclidean Distance between two points.
            D = pdist([pix_map(2) pix_map(1); tmp_orig(1) tmp_orig(2)],'euclidean');
            
            % If distance is less than p then, it is inliner.
            if (D < p)
                inliners = [inliners; i];
            end
        end
    end
 
    % It is the best estimate, if it has the most number of inliners.
    if (best_fit_est(1) < size(inliners,1))
        selected_inliners = inliners;
        best_fit_est = [size(inliners,1) ransac_pts(1) ransac_pts(2) ransac_pts(3)];
    end
end


%% PART. 13
close all;
right_t = [zeros(175,size(right,2)); right];
concat = [left right_t]; % concat images

% find coordinates
for i=1:size(selected_inliners,1)
    p1 = f1(:,filtered_match(selected_inliners,1));
    p2 = f2(:,filtered_match(selected_inliners,2));
end
p2(1,:) = p2(1,:) + size(left,2);
p2(2,:) = p2(2,:) + 175;

% draw!
figure, imshow(concat);
h1 = vl_plotframe(p1);
h2 = vl_plotframe(p2);
set(h1,'color','k','linewidth',3);
set(h1,'color','y','linewidth',2);
hold on;
for i = 1 : size(selected_inliners,1)
    plot([p1(1,i),p2(1,i)],[p1(2,i),p2(2,i)],'Color','r','LineWidth',1);
end
pause;

%% PART. 14
% Using all the inliers of the best transformation found using RANSAC (i.e., the one with the most inliers),
% compute the ?nal transformation with least-squares.

% Now we have the best estmate points and we can grap the transformation.
ransac_pts = [best_fit_est(2) best_fit_est(3) best_fit_est(4)];

A = [];
b = [];

for i = 1: size(selected_inliners,1)
    % Get the coordinates.
    p1f1 = f1(1:2, filtered_match(selected_inliners(i),1));
    p1f2 = f2(1:2, filtered_match(selected_inliners(i),2));
    
    % Estimate the affine transformation
    A = [A; ...
        p1f1(2),p1f1(1),0,0,1,0; ...
        0,0,p1f1(2),p1f1(1),0,1];
    
    b = [b; p1f2(2); p1f2(1)];
end

x = A\b;

% Now, we have our unknowns.
T = [x(1) x(2); x(3) x(4)];
c = [x(5); x(6)];

%% PART. 15
% Affine transformation.
T = maketform('affine', [x(1), x(2), 0; x(3), x(4), 0; x(5), x(6), 1]);
left_t = imtransform(left, T);
intersection = 1160; % manual padding for better looking.

composition = [ left_t zeros(size(left_t,1),intersection) ]; % put left ONTO right.
right_t = [ zeros(200, size(right,2)); right ]; 
right_window = size(right,2) - intersection; % intersecting window size.
left_window = size(composition,2) - intersection; % again.
intersect_pt = left_window-right_window-1; % Where it intersects.
composition(1:2600, left_window:end) = ...
    right_t(1:2600,right_window:end); % ONTO
composition(1:2600,intersect_pt:left_window) = ...
    composition(1:2600,intersect_pt:left_window)/2 + ...
    right_t(1:2600,1:right_window+2)/2;

close all;
figure, imshow(composition), pause;

%% PART. 16
orig_left = imread('building_left.jpg');
orig_right = imread('building_right.jpg');
left = im2single(rgb2gray(orig_left));
right = im2single(rgb2gray(orig_right));

[f1,d1] = vl_sift(left);
[f2,d2] = vl_sift(right);
[matches, scores] = vl_ubcmatch(d1, d2);

filtered_match = zeros(size(matches,2), 2);
index = 1;
for i = 1:size(matches,2) 
    % These images don't have much descriptors; i.e. about 300. 
    % So I decided to use them all to try.
    filtered_match(index,1) = matches(1,i);
    filtered_match(index,2) = matches(2,i);
    index = index + 1;
end
filtered_match = filtered_match(1:index-1,:);
best_fit_est = [0 0 0 0];
selected_inliners = [];

for N = 1:100
    ransac_pts = randperm(size(filtered_match,1),3);
    p1f1 = f1(1:2, filtered_match(ransac_pts(1),1));
    p1f2 = f2(1:2, filtered_match(ransac_pts(1),2));
    p2f1 = f1(1:2, filtered_match(ransac_pts(2),1));
    p2f2 = f2(1:2, filtered_match(ransac_pts(2),2));
    p3f1 = f1(1:2, filtered_match(ransac_pts(3),1));
    p3f2 = f2(1:2, filtered_match(ransac_pts(3),2));
    A = [p1f1(2),p1f1(1),0,0,1,0; ...
        0,0,p1f1(2),p1f1(1),0,1; ...
        p2f1(2),p2f1(1),0,0,1,0; ...
        0,0,p2f1(2),p2f1(1),0,1; ...
        p3f1(2),p3f1(1),0,0,1,0; ...
        0,0,p3f1(2),p3f1(1),0,1];
    b = [p1f2(2); p1f2(1); p2f2(2); p2f2(1); p3f2(2); p3f2(1)];
    x = A\b;
    T = [x(1) x(2); x(3) x(4)];
    c = [x(5); x(6)];
    
    p = 10; % To Keep the selected inliners in around 100-300 range.
    inliners = [];
    for i = 1 : size(filtered_match,1)
        if (ismember(i,ransac_pts) == 0)
            tmp_pix = f1(1:2, filtered_match(i,1));
            pix = [tmp_pix(2); tmp_pix(1)];
            pix_map = T*pix+c;
            
            tmp_orig = f2(1:2, filtered_match(i,2));
            
            D = pdist([pix_map(2) pix_map(1); tmp_orig(1) tmp_orig(2)],'euclidean');
            
            % If distance is less than p then, it is inliner.
            if (D < p)
                inliners = [inliners; i];
            end
        end
    end
 
    % It is the best estimate, if it has the most number of inliners.
    if (best_fit_est(1) < size(inliners,1))
        selected_inliners = inliners;
        best_fit_est = [size(inliners,1) ransac_pts(1) ransac_pts(2) ransac_pts(3)];
    end
end
ransac_pts = [best_fit_est(2) best_fit_est(3) best_fit_est(4)];
A = [];
b = [];
for i = 1: size(selected_inliners,1)
    p1f1 = f1(1:2, filtered_match(selected_inliners(i),1));
    p1f2 = f2(1:2, filtered_match(selected_inliners(i),2));
    A = [A; ...
        p1f1(2),p1f1(1),0,0,1,0; ...
        0,0,p1f1(2),p1f1(1),0,1];
    b = [b; p1f2(2); p1f2(1)];
end
x = A\b;
T = maketform('affine', [x(1), x(2), 0; x(3), x(4), 0; x(5), x(6), 1]);
left_t = imtransform(left, T);
intersection = 490; % manual padding for better looking.

composition = [ left_t zeros(size(left_t,1),intersection) ];
right_window = size(right,2) - intersection;
left_window = size(composition,2) - intersection;
intersect_pt = left_window-right_window-1;
composition(1:600, left_window:end) = ...
    right(1:600,right_window:end);
composition(1:600,intersect_pt:left_window) = ...
    composition(1:600,intersect_pt:left_window)/2 + ...
    right(1:600,1:right_window+2)/2;

close all;
figure, imshow(composition);







