% Copyright (c) 2008, Chris Marsh
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the University of Saskatchewan nor the names 
%       of its contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

clear all
close all

%basin or sun
viewpoint='sun';

load square_nodes_5m.csv

x=square_nodes_5m(:,1);
y=square_nodes_5m(:,2);
z=square_nodes_5m(:,3);

load boundary_nodes.csv
xb=boundary_nodes(:,1);
yb=boundary_nodes(:,2);
zb=boundary_nodes(:,3);

num_nodes = numel(x);


proj_x=zeros(num_nodes,1);
proj_y=zeros(num_nodes,1);
proj_z=zeros(num_nodes,1);

lat = 50.960873;
long = -115.187890  ;
alt = 0.0;


% datestr(time,'yyyy/mm/dd HH:MM:SS')
%utc offset
UTC_offset = 7;
tstart = datenum('2010/09/15 6:00:00');
tend = datenum('2010/09/15 18:00:00');
tstart=datestr(addtodate(tstart,UTC_offset,'hour'), 'yyyy/mm/dd HH:MM:SS');
tend=datestr(addtodate(tend,UTC_offset,'hour'), 'yyyy/mm/dd HH:MM:SS');

%correct for time zone
% UTC=datestr(addtodate(time,UTC_offset,'hour'), 'yyyy/mm/dd HH:MM:SS');
step_size = 30;
step_type = 'minute';
    
time = tstart;

frame=1;

tri=DelaunayTri(x,y);
color=zeros(size(tri.Triangulation));
hs=[];

%spatial grid
%number of cells in x direction
sg_nx=3;

if strcmp(viewpoint,'basin')
    hp=plot3(xb, yb, zb,'o','Color','black');
    %do this to force the plot to go to monitor 1 in a dual monitor setup. This
    %is needed for getframe, which won't work on a 2nd monitor. Seriously.
    movegui(hp)

    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
    axis tight
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    campos(  1.0e+006 .*[ 0.6651    5.6380    0.0080] )
    hold on
else
    %for "as the sun"
    ff=figure;
    movegui(ff)
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
end

%triangle centers
% triC=zeros(length(tri.Triangulation),3);
% shadows=zeros(length(tri.Triangulation),1);
% 
% for i=1:length(tri.Triangulation)
%      triC(i,1)=1/3*(x(tri.Triangulation(i,1))+x(tri.Triangulation(i,1))+x(tri.Triangulation(i,1)));
%      triC(i,2)=1/3*(y(tri.Triangulation(i,2))+y(tri.Triangulation(i,2))+y(tri.Triangulation(i,2)));
%      triC(i,3)=1/3*(z(tri.Triangulation(i,3))+z(tri.Triangulation(i,3))+z(tri.Triangulation(i,3)));
% end

    
while datenum(time, 'yyyy/mm/dd HH:MM:SS') <= datenum(tend, 'yyyy/mm/dd HH:MM:SS')
    
    [Az El] = SolarAzEl(time,lat,long,alt);

    %convert to rads
    A  = Az * pi/180;
    E  = El * pi/180;

    %eqns (6) & (7) in Montero
    z0 = pi-A;
    q0 = pi/2 - E;

    %projection
    K = [cos(z0) sin(z0) 0; -cos(q0)*sin(z0) cos(q0)*cos(z0) sin(q0);sin(q0)*sin(z0) -cos(z0)*sin(q0) cos(q0)];

    for i = 1:num_nodes

        coord = [x(i); y(i); z(i)];
        coord = K*coord;
        proj_x(i) = coord(1);
        proj_y(i) = coord(2);
        proj_z(i) = coord(3);%should this be -coord(3) ?? -z axis ?? no...
    end
    
    %OBB computation 
%     -----------------------
    
%         4       top      3
%         +--------+--------+
%         |        |        |
%         |        |        |
%   left  +--------+--------+    right
%         |        |        |
%         |        |        |
%         +-----------------+
%         1     bottom      2
% +y
% ^
% |
% +-> +x

    %only need 1:4 of these 5x1 vectors bbx, bby
    [bbx,bby,area,per]=minboundrect(proj_x,proj_y);
    tbb=[bbx(1:4),bby(1:4)]; %ignore the last point because it is out of order. Indecies follow the above naming convention
    
  
    m=(bby(2)-bby(1))/(bbx(2)-bbx(1));
    step=(bbx(3)-bbx(4))/sg_nx;
    for i=1:sg_nx
        xpos.bottom=bbx(1)+step*i;
        xpos.top=bbx(4)+step*i;
        
        x_mp.bottom(i).x=xpos.bottom;
        x_mp.bottom(i).y=m*(xpos.bottom-bbx(1))+bby(1); %2pt line eqn
        
        
        x_mp.top(i).x=xpos.top;
        x_mp.top(i).y=m*(xpos.top-bbx(4))+bby(4); %2pt line eqn
        
        plot(x_mp.bottom(i).x,x_mp.bottom(i).y,'o','color','blue')
        plot(x_mp.top(i).x,x_mp.top(i).y,'o','color','blue')
    end
     

    
    
    rectangles=cell(sg_nx,1);
    for i=1:sg_nx
        %first rect
       if i ==1
           rectangles{1}=[[bbx(1),bby(1)];...%bottom left
                         [x_mp.bottom(1).x,x_mp.bottom(1).y];...%bottom right mid point
                         [x_mp.top(1).x,x_mp.top(1).y];... %top right mid point
                         [bbx(4),bby(4)];...%top right
                         [bbx(1),bby(1)]]; 
%        last rect
       else if i==sg_nx
           rectangles{i}= [[x_mp.bottom(i-1).x,x_mp.bottom(i-1).y];...
                           [bbx(2),bby(2)];...
                           [bbx(3),bby(3)];...
                           [x_mp.top(i-1).x,x_mp.top(i-1).y];...
                           [x_mp.bottom(i-1).x,x_mp.bottom(i-1).y]];
           else
             rectangles{i}=[[x_mp.bottom(i-1).x,x_mp.bottom(i-1).y];...
                           [x_mp.bottom(i).x,x_mp.bottom(i).y];...
                           [x_mp.top(i).x,x_mp.top(i).y];...
                           [x_mp.top(i-1).x,x_mp.top(i-1).y];...
                           [x_mp.bottom(i-1).x,x_mp.bottom(i-1).y]];
           end
       end
           
    end
    
    
%     new_tri=[tri.Triangulation(:,:) triC(:,3)];
%     new_tri=sort(new_tri,4,'descend');
%     fprintf('building shadows');
%     parfor i=1:length(tri.Triangulation)
%         for j=i:length(tri.Triangulation)
%             %  z' of our ith triangle is behind our jth triangle
%             if triC(i,3) < triC(j,3)
%                if inside_triangle(triC(i,1:2), ...
%                        [proj_x(tri.Triangulation(j,1)) proj_y(tri.Triangulation(j,1))],...
%                        [proj_x(tri.Triangulation(j,2)) proj_y(tri.Triangulation(j,2))],...
%                        [proj_x(tri.Triangulation(j,3)) proj_y(tri.Triangulation(j,3))]) == 1 %true
%                   shadows(i)=1; %ith triangle is in shadow 
%                end
%             end
%         end
%         
%     end

% length(tri.Triangulation)
%     for i=1:length(tri.Triangulation)
%         for j=1:length(tri.Triangulation)
%             %  z' of our ith triangle is behind our jth triangle
%             if triC(i,3) < triC(j,3)
%                if inside_triangle(triC(i,1:2), ...
%                        [proj_x(tri.Triangulation(j,1)) proj_y(tri.Triangulation(j,1))],...
%                        [proj_x(tri.Triangulation(j,2)) proj_y(tri.Triangulation(j,2))],...
%                        [proj_x(tri.Triangulation(j,3)) proj_y(tri.Triangulation(j,3))]) == 1 %true
%                   shadows(i)=1; %ith triangle is in shadow 
%                end
%             end
%             j
%             
%         end 
%         
%     end
    r=1;
    sun = [r*cos(E)*cos(A);r*cos(E)*sin(A);r*sin(E)];

    if exist('ht','var')==1
        delete(ht.th);
    end

    
%testing code that plots the rectangles
% ------------------------------------------------'
% Plots the sub rect points
    % for i=1:4
    %     hold on
    %    plot(tbb(i,1),tbb(i,2),'o');text(tbb(i,1)+5,tbb(i,2)+5,num2str(i)); 
    %    pause
    % end

    %plot bounding rectangle
%     plot(bbx(:),bby(:),'color','red','linewidth',5);
% 
% Plots the sub rects
%     hold on
%         for i=1:sg_nx
%            plot(rectangles{i}(:,1),rectangles{i}(:,2),'color','blue','linewidth',5);
%         end
%     hold off
%     pause
        
%----------------------------------------
    if  exist('p','var')==0
        if strcmp(viewpoint,'basin')
            p = patch( ...
            'Vertices',[x(:) y(:) z(:)], ...
            'Faces', tri.Triangulation, ...
            'facevertexcdata',proj_z(:),...
            'facecolor','flat',...
            'edgecolor','black');
        else
    
%         for "as the sun"
            p = patch( ...
                'Vertices',[proj_x(:) proj_y(:)], ...
                'Faces', tri.Triangulation, ...
                'facevertexcdata',proj_x(:),...
                'facecolor','flat',...
                'edgecolor','black');
        end

    else
        if strcmp(viewpoint,'basin')
             set(p,'facevertexcdata',proj_z(:));
        else
            %for "as the sun"
            set(p,'Vertices',[proj_x(:) proj_y(:) ],'facevertexcdata',proj_x(:));
  
            axis tight

        end
        refreshdata
    end
    ht= mtit(datestr(addtodate(datenum(time, 'yyyy/mm/dd HH:MM:SS'),-UTC_offset,'hour'), 'yyyy/mm/dd HH:MM:SS'), 'fontsize',14) ;
  

%     F(frame)=getframe(gcf);
%     export_fig( num2str(frame),'-png');
    frame=frame+1;
     
    
    
    time = datestr(addtodate(datenum(time,'yyyy/mm/dd HH:MM:SS'),step_size,step_type),'yyyy/mm/dd HH:MM:SS');
   pause
end
hold off


