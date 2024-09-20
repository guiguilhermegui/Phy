 %%% PARTISAN (PARTicle Shape ANalyzer) v2.0 -- Parameterizing volcanic rock fragments
%   after Dellino and La Volpe (1996); Durig et al. (2012); Cioni et al. (2014);
%         Leibrandt and Le Pennec (2015); Liu et al. (2015); and Schmith et al. (2017).
%
%  Copyright (C) 2021 by Tobias Durig and M. Hamish Bowman
%
%  This is an updated version of PARTISAN v1.0 by M. Hamish Bowman and Tobias Durig.
%
%  Please cite Durig, T. et al. (2018) "PARTIcle Shape ANalyzer PARTISAN -
%    an open source tool for multi-standard two-dimensional particle
%    morphometry analysis", An. Geophys., 61, 31. doi: 10.4401/ag-7865
%
% 
%  This program is free software licensed under the GPL (>=v3).
%  Read the GPS-3.TXT file that comes with this software for details.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This program is free software; you can redistribute it and/or modify it
%  under the terms of the GNU General Public License as published by the
%  Free Software Foundation; either version 3 of the License, or (at your
%  option) any later version.
%  
%  Parts of PARTISAN are not copyright by the PARTISAN development team.
%  The original authors hold the copyrights and you have to abide to their
%  licensing terms where noted. See the headers of the respective .m scripts
%  for details.
%  
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License (GPL) for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, see <http://www.gnu.org/licenses/>.
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input is in the form of a series of silhouette image files of any format
%       known to MATLAB's image reading function.
%  Output is to a flat ASCII file named 'particle_stats_output.csv', as parameters
%       displayed in the command window, and (optionally) to a series of figures
%       displayed on the screen.
%
%  Requires MATLAB's Imagery toolbox.  Tested with MATLAB version R2016b.
%    (GNU Octave support is TODO: bwconvhull() is missing)
%
%  Adjust the 'do_plots' and 'img_files' below as required.
%  The 'img_files' may contain wildcards and at your option the full or relative
%    path. e.g., '*.tif', or 'example1/*.bmp' or 'C:\Image Files\*.png'.
%
%  Cell-center vs. grid-confluence conventions are not always consistent between e.g.
%   area and perimeter support functions, but we do our best in an imperfect world.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
do_plots = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cd 'U:\PHY\PARTISAN'
cd U:\PHY\PARTISAN\Examples
%cd 'U:\PHY\PARTISAN\Examples\RDI_ArcelorMittal_SEM\RDI'
img_files = dir('*.tif');      %%%%%%%%%%%  <-- adjust as needed  %%%%%%%%%%%%
%img_files  dir('*.jpg');
%img_files = [dir('*.bmp')  dir('*.tif')  dir('*.TIF')  dir('*.png')  dir('*.jpg')];
%img_files = [dir('*.jfif')  dir('*.png')  dir('*.jpg')];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fnames = {img_files.name};
if(isfield(img_files, 'folder'))  % doesn't exist in older versions of Matlab
   fpath = {img_files.folder};
   fpath = fpath{1};  %%% se for um arquivo só nas linhas 60, pode comentar
else
   fpath = pwd;   % better than nothing :-/
end
output_is_open = false;


for j=1:length(fnames)
   filename = fnames{j};

   data = imread([fpath '/' filename]);
   data=data(:,:,1); %% Prof Giuliano   tirar mais de uma camada
   %data = imbinarize(data) %%%% teste guilherme 16/06 15h
   if(~output_is_open)
      %% prepare output file
      fid = fopen('particle_stats_output.csv', 'w');
      if(fid < 0)
         disp('ERROR: could not open the output file <particle_stats_output.csv>.')
         break
      else
         output_is_open = true;
      end
      
      fprintf(fid, ['PARTISAN (PARTicle Shape ANalyzer) particle parameters -- run ' datestr(now) '\n\n']);

      fprintf(fid, ',Basic metrics,,,,,,,,,,,,,,,,,,,,IPA: Dellino and La Volpe (1996),,,,,');
      fprintf(fid, 'Cioni et al (2014),,,,,Leibrandt and Le Pennec (2015),,,,,,,');
      fprintf(fid, 'Liu et al (2015),,,,,Schmith et al (2017),,,,\n');   
      fprintf(fid, ['filename,perimeter,particle area,width,breadth,perimeter of circle with area A,' ...
	      'maximum intercept,mean intercept perpendicular,' ...
	      'max. proj. length,orthogonal proj. length,Heywood diameter,convex hull perimeter,' ...
	      'convex hull area,min ellipse perimeter,major best fit ellipse,minor best fit ellipse,' ...
	      'diameter bounding circle,minimum Feret,orthogonal Feret,maximum Feret,,' ...
	      'circularity,rectangularity,compactness,elongation,,' ...
	      'circularity,ellipse aspect ratio,convexity,solidity,,' ...
	      'circularity,elongation, aspect ratio,CE diameter,convexity,solidity,,' ...
	      'form factor,axial ratio,convexity,solidity,,' ...
	      'circularity,rectangularity,form factor,Feret aspect ratio,Feret aspect ratio,regularity parameter\n']);
      fprintf(fid, [',p,A,w,b,c,a,m,L_b,W_b,d_H,p_CP,A_CP,e_CE,L_maj,L_min,d_BC,' ...
	      'l_F,w_F,d_F,,Circ_DL,Rec_DL,Com_DL,Elo_DL,,Circ_CI,AR_CI,Con_CI,Sol_CI,,' ...
	      'Circ_LL,Elo_LL,AR_LL,d_H,Con_LL,Sol_LL,,FF,AR_LI,Con_LI,Sol_LI,,' ...
	      'Circ_SC,Rec_SC,FF,AR_F,AR_SC,Reg\n']);
   end


   %%% test sim data:  749x949
   if(false)
     test_blank = zeros(749,949);

     test_3x3 = test_blank;
     test_3x3(350:352, 500:502) = ones(3,3) * 255;

     test_2x4 = test_blank;
     test_2x4(350:351, 500:503) = ones(2,4) * 255;

     %data = test_3x3;
     data = test_2x4;
   end
   %%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BASIC STATS ON RAW IMAGE

im_dims = size(data);
tot_num_cells = im_dims(1) * im_dims(2);
cell_size = 1;   % apply scaling and units as needed

if( max(max(data) == 1) )
  dthresh = 0.5;
elseif(max(max(data) == 255) )
  dthresh = 128;
else
  dthresh = max(max(data))./2; %villa mod -> OLD original  dthresh = 128;
  disp([      'FIXME: Max value of the data is non-standard (' ...
       num2str(max(max(data))) '). Using ' num2str(dthresh) ...
       ' as the threshold between presence and absence.' ...
       'Último aviso Villa check -> apenas pra dizer que a cor máxima não é branco (255) '])
end


if(data(1,1) > dthresh)
   disp(['***' filename ' data is inverted! inverting back..']) 
   if( max(max(data) == 1) )
      data = 1 - data;
   elseif(max(max(data) == 255) )
       data = 255 - data;
   else
       data = 255 - data;
   end
end

num_cells = length(find(data > dthresh));   % num cells with data
particle_area = bwarea(data);

A = num_cells * cell_size;

% Matlab's image processing toolbox's version of things
%  https://au.mathworks.com/help/images/ref/regionprops.html
data_bool = data > dthresh;
imgstats = regionprops(data_bool, 'Area', 'BoundingBox', 'Centroid', ...
  'ConvexArea', 'ConvexHull', 'ConvexImage', 'Eccentricity', ...
  'EquivDiameter', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', ...
  'Perimeter', 'Solidity');


if(length(imgstats) > 1)
   % if there is more than one area (island), only consider the largest of them.
   maxArea = 0;
   biggestArea = NaN;
   for i = 1:length(imgstats)
      if (imgstats(i).Area > maxArea)
         biggestArea = i;
	 maxArea = imgstats(i).Area;
      end
   end
   imgstats = imgstats(biggestArea);
   clear maxArea biggestArea
end


%% find convex hull of the image
khull = bwconvhull(data);
% combine the null background with convex hull and positive image data (for plotting)
kmix = zeros(size(data));
kmix(find(khull)) = 127;
kmix(find(data > dthresh)) = 255;


% find area of the convex hull image cells
Ahull = length(find(khull)) * cell_size;


% find the minimum bounding box
[Ys, Xs] = find(data);
cc = minBoundingBox([Xs Ys]');  % corner_coords


% find the perimeter (fractal problem; angularized centroids, and other caveats apply)
Bdy = bwboundaries(data);


if(do_plots)
   figure
   %set(gcf, 'Position', [30 50 1500 450])   % for 3 plots in a row
   set(gcf, 'Position', [10 40 1000 800])  % for 2x2 plots

   subplot(2, 2, 1)
   %title('Raw image', 'FontWeight', 'normal')
   %imagesc(data)
   imagesc(kmix)
   colormap([1 1 1; 0 0 1; 0.5 0.5 0.5])%!!!
   axis equal
   axis tight
   grid on
   hold on
   plot(cc(1,[1:end 1]), cc(2,[1:end 1]), 'k','LineWidth', 1.5) %!!!
   % debug:
   %for i=1:4
   %   text(cc(1,i), cc(2,i), num2str(i))
   %end

   % perimeter
   for k = 1:length(Bdy)
      if(k==1)
         bColor = 'y';
      else
         bColor = 'k';
      end
      boundary = Bdy{k};
      plot(boundary(:,2), boundary(:,1), bColor, 'LineWidth', 2)
   end
end


d = diff([Bdy{1}(:,1) Bdy{1}(:,2)]);
perimeter = sum(sqrt(sum(d.^2, 2)));
p = perimeter;


% find the perimeter of the convex hull
BdyHull = bwboundaries(khull);
dHull = diff([BdyHull{1}(:,1) BdyHull{1}(:,2)]);
pHull = sum(sqrt(sum(dHull.^2, 2)));
clear dHull


b1 = sqrt((cc(1, 2) - cc(1, 1))^2 + (cc(2, 2) - cc(2, 1))^2 );    % breadth
%same b2 = sqrt((cc(1, 4) - cc(1, 3))^2 + (cc(2, 4) - cc(2, 3))^2 )    % breadth

w1 = sqrt((cc(1, 3) - cc(1, 2))^2 + (cc(2, 3) - cc(2, 2))^2 );   % width
%same w2 = sqrt((cc(1, 4) - cc(1, 1))^2 + (cc(2, 4) - cc(2, 1))^2 )   % width


% check that breadth is the longer of the two sides, if not reverse them
%  and rotate things by 90 degrees
if(b1 > w1)
  b = b1;
  w = w1;
  ExtraRot = 0;
else
  b = w1;
  w = b1;
  ExtraRot = pi/2;
end


feret_minor = w;
w_f = b;


% compactness = A / (b*w);

% A = pi * r^2
r = sqrt(A/pi);
c = pi * 2*r;   % perimeter of a circle with area A

% Heywood diameter
dH = r * 2;

%% find and plot the center of gravity
cog_xy = cog_coord(data, dthresh);

%%% rotate box so that longest axis is horizontal
%dx = cc(1,4) - cc(1,3);
dy = cc(2,4) - cc(2,3);

%sin theta = opp/hyp = dy/b
theta = real(asin(dy/b1)) + ExtraRot;
%theta = acos(dx/b);


%% minimum bounding circle
[minCirc.Ys, minCirc.Xs] = find(data);
[minCirc.center, minCirc.radius] = minBoundCircle(minCirc.Xs, minCirc.Ys, false);
minCirc.diameter = minCirc.radius * 2;
d_BC = minCirc.diameter;


if(do_plots)
   % plot the equiv area circle
   viscircles([mean(cc(1,:)) mean(cc(2,:))], r, ...
      'Linewidth', 2, 'color', [0.9 0.4 0.1]);%!!!

   % and bbox's geometric center
   plot(mean(cc(1,:)), mean(cc(2,:)), '+', 'color', [0.9 0.4 0.1],'LineWidth', 1.5)%!!!

   plot(cog_xy(1), cog_xy(2), 'k+')

   viscircles(minCirc.center, minCirc.radius, 'Linewidth', 1.5, 'color', [.4 .8 .4]);%!!!

   text(cc(1, 3), cc(2,3), [' ' num2str(abs(theta)*180/pi, '%.1f') char(176)], ...
       'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
       'fontsize', 9, 'rotation', -1*theta*180/pi)

   text(minCirc.center(1) + minCirc.radius * cos(pi/4), ...
        minCirc.center(2) - minCirc.radius * sin(pi/4), '\itBC', ...
	'VerticalAlign', 'bottom', 'color', [.4 .8 .4])%!!!

   text(mean(cc(1,:)) + r * cos(pi/4), ...
        mean(cc(2,:)) - r * sin(pi/4), '\itc', ...
	'VerticalAlign', 'bottom', 'color', [0.9 0.4 0.1])%!!!

   if(ExtraRot)
      text((cc(1,3) + cc(1,4))/2, (cc(2,3) + cc(2,4))/2, '\itw', ...
           'color', 'k', 'VerticalAlign', 'top')%!!!
      text((cc(1,2) + cc(1,3))/2, (cc(2,2) + cc(2,3))/2, '\itb ', ...
           'color', 'k', 'HorizontalAlign', 'right', 'VerticalAlign', 'top')
   else
      text((cc(1,3) + cc(1,4))/2, (cc(2,3) + cc(2,4))/2, ' \itb', ...
           'color', 'k', 'VerticalAlign', 'top', 'HorizontalAlign', 'left')%!!!
      text((cc(1,2) + cc(1,3))/2, (cc(2,2) + cc(2,3))/2, '\itw ', ...
           'color', 'k', 'HorizontalAlign', 'right', 'VerticalAlign', 'top')%!!!
   end

   %% filename as figure title
   % three plots in a row: try placement at 1.775, 1.08
   text(1.3, 1.15, filename, 'Interpreter', 'none', 'FontSize', 12, ...
       'HorizontalAlignment', 'center', 'Units','normalized')
   set(gcf, 'Name', filename)   % for figure window's title bar
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ROTATED TO ALIGN WITH MINIMUM BOUNDING BOX

% rotated so that bbox is "landscape" orientation
data_rot_bbox = imrotate(data, theta * 180/pi, 'bilinear');
[Ys, Xs] = find(data_rot_bbox);
cc_rot = minBoundingBox([Xs Ys]');
crop_x0 = round(min(cc_rot(1,:)));
crop_y0 = round(min(cc_rot(2,:)));
crop_w = round(max(cc_rot(1,:)) - min(cc_rot(1,:)));
crop_h = round(max(cc_rot(2,:)) - min(cc_rot(2,:)));
data_rot_crop_bbox = imcrop(data_rot_bbox, [crop_x0 crop_y0 crop_w crop_h]);


%%%% determine "a": max segment length parallel to bbox's longer dimension
max_segment.value = 0;
max_segment.row = NaN;
max_segment.part = NaN;
max_segment.start = 0;

% scan through row by row and part by part
for i = 1:size(data_rot_crop_bbox, 1)
   % skip rows with just a shade of sample in it
   hits = find(data_rot_crop_bbox(i,:) > dthresh);
   if(length(hits) > 0)
      jumps = find(diff(hits) > 1);  % index within hits where the next record jumps
      if(length(jumps) == 0)
         len = 1+ max(hits) - min(hits);
         if(len > max_segment.value)
	    max_segment.value = len;
	    max_segment.row = i;
	    max_segment.part = 1;
	    max_segment.start = min(hits);
          end
      else
         for j=1:length(jumps)
            if(j==1)
	       start_seg = min(hits);
	    else
	       start_seg = hits(jumps(j-1)+1);
            end

            len = 1+ hits(jumps(j)) - start_seg;

            if(len > max_segment.value)
	       max_segment.value = len;
	       max_segment.row = i;
	       max_segment.part = j;
	       max_segment.start = start_seg;
            end
	 end
      end
   end
end

a = max_segment.value;

%%% m = mean intercept perpendicular to "a"
chord_sum = 0;
num_chords = 0;
% scan through column by column
for i = 1:size(data_rot_crop_bbox, 2)
   % skip columns with just a shade of sample in it
   hits = find(data_rot_crop_bbox(:,i) > dthresh);
   if(length(hits) > 0)
      len = 1 + max(hits) - min(hits);
      chord_sum = chord_sum + len;
      num_chords = num_chords + 1;
   end
end

% chord_sum - A    % (should be about the same)
m = chord_sum / num_chords;

%% find the center of gravity/mass/barycenter
cog_xyBBox = cog_coord(data_rot_crop_bbox, dthresh);


% plotting
if(do_plots)
   subplot(222)
   %title('Bounding box aligned')  % commented out here; added at the end
   imagesc(data_rot_crop_bbox)
   axis equal
   axis tight
   grid on
   hold on

   plot([max_segment.start max_segment.start+max_segment.value], ...
        [max_segment.row max_segment.row], 'b','linewidth',1.5)
   plot([mean(xlim) mean(xlim)], [mean(ylim)+m/2 mean(ylim)-m/2], 'b','LineWidth',1.5)
   % bbox's geometric center
   hCoBB = plot(size(data_rot_crop_bbox,2)/2, size(data_rot_crop_bbox,1)/2, ...
               '+', 'color', [0.9 0.4 0.1],'LineWidth',1.5);%!!!
   hCoG = plot(cog_xyBBox(1), cog_xyBBox(2), 'k+');


   if(max(ylim) - max_segment.row < 25)
      tmp.valign = 'bottom';
   else
      tmp.valign = 'top';
   end
   text((max_segment.start+max_segment.value)/5 + max_segment.start, ...
        max_segment.row, '\ita', 'VerticalAlign', tmp.valign)
   clear tmp

   text(mean(xlim), mean(ylim) + m/4, ' \itm', 'HorizontalAlign', 'left')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ROTATED TO ALIGN WITH FERET DIAMETER

% rotate in 1/2 degree steps. start with theta0 and look a bit more than 45 deg +/-
max_intercept.value = 0;
max_intercept.row = NaN;
max_intercept.theta = NaN;


% speed things up for large images  (perhaps not good to alter method between particles?)
if(prod(size(data)) > 512^2)
   stepsize = 1.0;
else
   stepsize = 0.5;
end

for adjust = -50 : stepsize : +50
   %disp((theta0 * 180/pi) + adjust)
   data_rot_test = imrotate(data, (theta * 180/pi) + adjust, 'bilinear');

   % scan through row by row
   for i = 1:size(data_rot_test, 1)
      % skip rows with just a shade of sample in it
      hits = find(data_rot_test(i,:) > dthresh);
      if(length(hits) > 0)
         len = 1+ max(hits) - min(hits);
         if(len > max_intercept.value)
            max_intercept.value = len;
            max_intercept.row = i;
	    max_intercept.theta = (theta * 180/pi) + adjust;
          end
      end
   end

end


data_rot = imrotate(data, max_intercept.theta, 'bilinear');
feret_rot_theta = max_intercept.theta;
clear data_rot_test max_intercept


[Ys, Xs] = find(data_rot);
cc_rot = minBoundingBox([Xs Ys]');
crop_x0 = round(min(cc_rot(1,:)));
crop_y0 = round(min(cc_rot(2,:)));
crop_w = round(max(cc_rot(1,:)) - min(cc_rot(1,:)));
crop_h = round(max(cc_rot(2,:)) - min(cc_rot(2,:)));

data_rot_crop = imcrop(data_rot, [crop_x0 crop_y0 crop_w crop_h]);
drc = data_rot_crop;



if(do_plots)
   subplot(223)
   %title('Feret aligned')

   imagesc(data_rot_crop)
   axis equal
   axis tight
   grid on
   hold on
   %plot(cc_rot(1,[1:end 1]), cc_rot(2,[1:end 1]),'r')
end


%%% maximum Feret distance = max horiz intercept
max_intercept.value = 0;
max_intercept.row = NaN;

% scan through row by row
for i = 1:size(drc, 1)
   % skip rows with just a shade of sample in it
   hits = find(drc(i,:) > dthresh);
   if(length(hits) > 0)
      len = 1+ max(hits) - min(hits);
      if(len > max_intercept.value)
         max_intercept.value = len;
         max_intercept.row = i;
       end
   end
end


feret_major = max_intercept.value;

%% find the center of gravity
cog_xyRot = cog_coord(drc, dthresh);


%% for plotting the bounding ellipse
ellip.a = imgstats.MajorAxisLength / 2;
ellip.b = imgstats.MinorAxisLength / 2;
% NOTE: this will often be close to the bbox's orientation, But It Is Not The Same.
ellip.theta = -1 * (imgstats.Orientation + feret_rot_theta) * pi/180;
ellip.x0 = cog_xyRot(1);
ellip.y0 = cog_xyRot(2);
ellip.Xs = ellip.x0 + ellip.a * cos(0 : 0.01 : 2*pi + 0.01);
ellip.Ys = ellip.y0 + ellip.b * sin(0 : 0.01 : 2*pi + 0.01);

rot.pts = [ellip.Xs; ellip.Ys];
rot.centered = repmat([ellip.x0; ellip.y0], 1, length(rot.pts));
rot.matrix = [cos(ellip.theta) -sin(ellip.theta); sin(ellip.theta) cos(ellip.theta)];
rot.pts_rot = rot.matrix * (rot.pts - rot.centered) + rot.centered;
rot.Xs = rot.pts_rot(1,:);
rot.Ys = rot.pts_rot(2,:);

rot.cardinal_Xs = nan(1,4);
rot.cardinal_Ys = nan(1,4);
rot.cardinal_Xs(1) = ellip.x0 + ellip.a * cos(0);
rot.cardinal_Xs(2) = ellip.x0 + ellip.b * cos(pi/2);
rot.cardinal_Xs(3) = ellip.x0 + ellip.a * cos(pi);
rot.cardinal_Xs(4) = ellip.x0 + ellip.b * cos(2*pi - pi/2);
rot.cardinal_Ys(1) = ellip.y0 + ellip.a * sin(0);
rot.cardinal_Ys(2) = ellip.y0 + ellip.b * sin(pi/2);
rot.cardinal_Ys(3) = ellip.y0 + ellip.a * sin(pi);
rot.cardinal_Ys(4) = ellip.y0 + ellip.b * sin(2*pi - pi/2);
rot.cardinal_pts = [rot.cardinal_Xs; rot.cardinal_Ys];
rot.card_cent = repmat([ellip.x0; ellip.y0], 1, length(rot.cardinal_pts));
rot.cardinal_pts_rot = rot.matrix * (rot.cardinal_pts - rot.card_cent) + rot.card_cent;
rot.cardinal_Xs = rot.cardinal_pts_rot(1,:);
rot.cardinal_Ys = rot.cardinal_pts_rot(2,:);


% containing (aka minimum bounding) ellipse and its perimeter
[mbe.C, mbe.coords] = minBoundEllipse(data_rot_crop, 0.005);
mbe.dist = diff([mbe.coords(1,:)' mbe.coords(2,:)']);
mbe.perimeter = sum(sqrt(sum(mbe.dist.^2, 2)));
ece = mbe.perimeter;
clear mbe.dist


if(do_plots)
   hits = find(drc(max_intercept.row,:) > dthresh);
   plot([min(hits) max(hits)], 0.5+[max_intercept.row max_intercept.row], 'b','LineWidth', 1.5)

   % and the containing ellipse:
   hCoBE = plot(mbe.C(1), mbe.C(2), '*', 'color', [1 0 1],'LineWidth', 1);
   plot(mbe.coords(1,:), mbe.coords(2,:), 'color', [1 0 1],'LineWidth', 2)

   % and the best-fit ellipse:
   plot(rot.Xs, rot.Ys, 'color', [.0 .9 .0], 'LineWidth', 2)

   plot([rot.cardinal_Xs(1) rot.cardinal_Xs(3)], ...
        [rot.cardinal_Ys(1) rot.cardinal_Ys(3)], ':', 'color', [.0 .9 .0],'LineWidth', 2)

   plot([rot.cardinal_Xs(2) rot.cardinal_Xs(4)], ...
        [rot.cardinal_Ys(2) rot.cardinal_Ys(4)], ':', 'color', [.0 .9 .0],'LineWidth', 2)

   plot(cog_xyRot(1), cog_xyRot(2), 'k+')


   text((cog_xyRot(1) + rot.cardinal_Xs(4))/2, ...
        (cog_xyRot(2) + rot.cardinal_Ys(4))/2, ' \itL_{min}', 'color', [.0 .9 .0])
   text((cog_xyRot(1) + rot.cardinal_Xs(1))/2, ...
        (cog_xyRot(2) + rot.cardinal_Ys(1))/2, ' \itL_{maj}', 'VerticalAlign', 'top', 'color', [.0 .9 .0])


   % contortions to make sure e_ce text always shows up on the top-right
   mbe.topYsIdx = find(mbe.coords(2,:) < mbe.C(2));
   pc90.x = max(mbe.coords(1,:)) * 0.9;
   [pc90.xmin, pc90.idx] = min(abs(mbe.coords(1, mbe.topYsIdx) - pc90.x));
   pc90.idx = mbe.topYsIdx(pc90.idx);

   text(mbe.coords(1, pc90.idx), mbe.coords(2, pc90.idx), ...
        '\ite_{ce}', 'color', [1 0 1], 'VerticalAlign', 'bottom')%[0.494 0.184 0.556]
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ROTATED TO ALIGN WITH BEST-FIT ELLIPSE

ellip.theta_raw = -1 * imgstats.Orientation;
data_rot_ellips = imrotate(data, ellip.theta_raw, 'bilinear');

[Ys, Xs] = find(data_rot_ellips);
crop_x0 = min(min(Xs));
crop_y0 = min(min(Ys));
crop_w = max(max(Xs)) - min(min(Xs));
crop_h = max(max(Ys)) - min(min(Ys));
data_rot_crop_ellips = imcrop(data_rot_ellips, [crop_x0 crop_y0 crop_w crop_h]);

cog_xyRotEllips = cog_coord(data_rot_crop_ellips, dthresh);

data_rot_crop_ellips_bool = data_rot_crop_ellips > dthresh;

Lb = max(find(max(data_rot_crop_ellips_bool))) - min(find(max(data_rot_crop_ellips_bool)));
Wb = max(find(max(data_rot_crop_ellips_bool'))) - min(find(max(data_rot_crop_ellips_bool')));
Lb_x0 = min(find(max(data_rot_crop_ellips_bool)));
Wb_y0 = min(find(max(data_rot_crop_ellips_bool')));


if(do_plots)
   subplot(224)

   imagesc(data_rot_crop_ellips)
   axis equal
   axis tight
   grid on
   hold on
   %plot(xlim, [cog_xyRotEllips(2) cog_xyRotEllips(2)], ':', 'color', [.32 .72 .32])
   %plot([cog_xyRotEllips(1) cog_xyRotEllips(1)], ylim, ':', 'color', [.32 .72 .32])

   plot([Lb_x0 Lb_x0+Lb], [cog_xyRotEllips(2) cog_xyRotEllips(2)], 'color', [.0 .9 .9], 'LineWidth', 1.5)
   plot([cog_xyRotEllips(1) cog_xyRotEllips(1)], [Wb_y0 Wb_y0+Wb], 'color', [.0 .9 .9], 'LineWidth', 1.5)
   plot(cog_xyRotEllips(1), cog_xyRotEllips(2), 'k+')

   text(cog_xyRotEllips(1) + Lb/6, cog_xyRotEllips(2), '\itL_b', 'VerticalAlign', 'top')
   text(cog_xyRotEllips(1), cog_xyRotEllips(2) + Wb/4, ' \itW_b', 'HorizontalAlign', 'left')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(do_plots)
   subplot(221)
   title('Raw image', 'FontWeight', 'normal')
   subplot(222)
   title('Bounding box aligned', 'FontWeight', 'normal')
   subplot(223)
   title('Feret aligned', 'FontWeight', 'normal')
   subplot(224)
   title('Best-fit ellipse aligned', 'FontWeight', 'normal')

   hL = legend([hCoBB hCoG hCoBE], 'Center of bounding box', 'Center of mass', ...
               'Center of bounding ellipse');
   set(hL, 'Position', [0.635 0.485 0.2 0.07])  %, 'Orientation', 'horizontal');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dellino and LaVolpe IPA (1996)
DL.Circ = p / c;
DL.Rec = p / (2*b+2*w);
DL.Com = A / (b*w);
DL.Elo = a / m;

% Cioni et al. (2014)
CI.Circ = (4*pi*A / p^2);
CI.AR = ellip.a / ellip.b;
CI.Con = ece / p;
CI.Sol = A / Ahull;

% Leibrandt and Le Pennec (2015)
LL.Circ = c / p;
LL.AR = Wb / Lb;
LL.Elo = 1 - LL.AR;
d_H = dH;
LL.Con = pHull / p;
LL.Sol = A / Ahull;

% Liu et al. (2015)
FF = (4*pi*A / p^2);
LI.AR= ellip.b / ellip.a;
LI.Con = pHull / p;
LI.Sol = A / Ahull;

% Schmith et al. (2017)
SC.Circ = 4*A / (pi * d_BC^2); 
SC.Rec = A / (b*w);
%FF = (4*pi*A / p^2);   % same as Liu et al.
SC.ARF =  feret_minor / w_f; % as defined in their Fig.1 
SC.AR = w_f / feret_minor; % displayed Feret aspect ratio (see their Fig.3)
Reg = SC.Circ * SC.Rec * FF;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elongation = a/m
% rectangularity = p/(2*w + 2*b)
% circularity = p/c

disp(['File: ' filename])
disp(' ')
disp('==== BASIC METRICS ===')
disp(['p: perimeter =                            ' num2str(p)])
disp(['A: particle area =                        ' num2str(A)])
disp(['w: width =                                ' num2str(w)])
disp(['b: breadth =                              ' num2str(b)])
disp(['c: perimeter of a circle with area A =    ' num2str(c)])
disp(['a: maximum intercept =                    ' num2str(a)])
disp(['m: mean intercept perpendicular =         ' num2str(m)])
disp(['L_b: max projected length       =         ' num2str(Lb)])
disp(['W_b: orthogonal projected length =        ' num2str(Wb)])
disp(['dH: Heywood diameter =                    ' num2str(dH)])
disp(['p_cp: convex hull perimeter =             ' num2str(pHull)])
disp(['A_cp: convex hull area =                  ' num2str(Ahull)])
disp(['e_ce: perimeter of min bounding ellipse = ' num2str(ece)])
disp(['L_maj: major axis of best fit ellipse =   ' num2str(2*ellip.a)])
disp(['L_min: minor axis of best fit ellipse =   ' num2str(2*ellip.b)])
disp(['d_BC: diameter of bounding circle =       ' num2str(d_BC)])
disp(['l_F: minimum Feret distance =             ' num2str(feret_minor)])
disp(['w_F: orthogonal Feret distance =          ' num2str(w_f)])
disp(['maximum Feret distance =                  ' num2str(feret_major)])
disp(' ')
disp('==== IPA: DELLINO and LA VOLPE (1996) ===')
disp(['Circ_DL: circularity =                    ' num2str(DL.Circ)])
disp(['Rec_DL: rectangularity =                  ' num2str(DL.Rec)])
disp(['Com_DL: compactness =                     ' num2str(DL.Com)])
disp(['Elo_DL: elongation =                      ' num2str(DL.Elo)])
disp(' ')
disp('==== CIONI et al. (2014) ===')
disp(['Circ_CI: circularity =                    ' num2str(CI.Circ)])
disp(['AR_CI: aspect ratio =                     ' num2str(CI.AR)])
disp(['Con_CI: convexity =                       ' num2str(CI.Con)])
disp(['Sol_CI: solidity =                        ' num2str(CI.Sol)])
disp(' ')
disp('==== LEIBRANDT and LE PENNEC (2015) ===')
disp(['Circ_LL: circularity =                    ' num2str(LL.Circ)])
disp(['Elo_LL: elongation =                      ' num2str(LL.Elo)])
disp(['AR_LL: aspect ratio =                     ' num2str(LL.AR)])
disp(['d_H: Heywood diameter =                   ' num2str(dH)])
disp(['Con_LL: convexity =                       ' num2str(LL.Con)])
disp(['Sol_LL: solidity =                        ' num2str(LL.Sol)])
disp(' ')
disp('==== LIU et al. (2015) ===')
disp(['FF: form factor =                         ' num2str(FF)])
disp(['AR_LI: axial ratio =                      ' num2str(LI.AR)])
disp(['Con_LI: convexity =                       ' num2str(LI.Con)])
disp(['Sol_LI: solidity =                        ' num2str(LI.Sol)])
disp(' ')
disp('==== SCHMITH et al. (2017) ===')
disp(['Circ_SC: circularity =                    ' num2str(SC.Circ)])
disp(['Rec_SC: rectangularity =                  ' num2str(SC.Rec)])
disp(['FF: form factor =                         ' num2str(FF)])
disp(['AR_F: Feret aspect ratio =                ' num2str(SC.ARF)])
disp(['AR_SC: Feret aspect ratio for plots =     ' num2str(SC.AR)])
disp(['Reg: regularity parameter =                ' num2str(Reg)])
disp(' ')
disp(' ')



fprintf(fid, [num2str(filename) ',' num2str(p) ',' num2str(A) ',' num2str(w) ',']);
fprintf(fid, [num2str(b) ',' num2str(c) ',' num2str(a) ',' num2str(m) ',']);
fprintf(fid, [num2str(Lb) ',' num2str(Wb) ',' num2str(dH) ',' num2str(pHull) ',']);
fprintf(fid, [num2str(Ahull) ',' num2str(ece) ',' num2str(2*ellip.a) ',']);
fprintf(fid, [num2str(2*ellip.b) ',' num2str(d_BC) ',' num2str(feret_minor) ',']);
fprintf(fid, [num2str(w_f) ',' num2str(feret_major) ',,' num2str(DL.Circ) ',' num2str(DL.Rec) ',']);
fprintf(fid, [num2str(DL.Com) ',' num2str(DL.Elo) ',,']);
fprintf(fid, [num2str(CI.Circ) ',' num2str(CI.AR) ',' num2str(CI.Con) ',']);
fprintf(fid, [num2str(CI.Sol) ',,' num2str(LL.Circ) ',' num2str(LL.Elo) ',']);
fprintf(fid, [num2str(LL.AR) ',' num2str(dH) ',' num2str(LL.Con) ',']);
fprintf(fid, [num2str(LL.Sol) ',,' num2str(FF) ',' num2str(LI.AR) ',']);
fprintf(fid, [num2str(LI.Con) ',' num2str(LI.Sol) ',,']);
fprintf(fid, [num2str(SC.Circ) ',' num2str(SC.Rec) ',' num2str(FF) ',']);
fprintf(fid, [num2str(SC.ARF) ',' num2str(SC.AR) ',' num2str(Reg)]);
fprintf(fid, '\n');

end %  filename loop


fclose(fid);


