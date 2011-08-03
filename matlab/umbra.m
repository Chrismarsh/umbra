% Copyright (c) 2011, Chris Marsh
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

function umbra()

%basin or sun
viewpoint='basin';

load square_nodes_2m.csv

x=square_nodes_2m(:,1); %#ok<NODEF>
y=square_nodes_2m(:,2);
z=square_nodes_2m(:,3);

load boundary_nodes.csv
xb=boundary_nodes(:,1); %#ok<NODEF>
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
tstart = datenum('2011/02/1 8:30:00');
tend = datenum('2011/02/1 19:00:00');
tstart=datestr(addtodate(tstart,UTC_offset,'hour'), 'yyyy/mm/dd HH:MM:SS');
tend=datestr(addtodate(tend,UTC_offset,'hour'), 'yyyy/mm/dd HH:MM:SS');

%correct for time zone
% UTC=datestr(addtodate(time,UTC_offset,'hour'), 'yyyy/mm/dd HH:MM:SS');
step_size = 15;
step_type = 'minute';
    
time = tstart;

frame=1;

%the triangulation
tri=DelaunayTri(x,y);
% Triangle_set=cell(tri.size,1);
    
%     
% color=zeros(size(tri.Triangulation));
% hs=[];

%spatial grid
%number of cells in x direction
sg_nx=3;
%pre allocation
rectangles=cell(sg_nx,1);

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

%triangluation
%x,y,z,proj_x,proj_y,proj_z,triCx,triCy,triCz,shadow
% Triangulation=zeros(max(tri.size),9);

shadows=zeros(length(tri.Triangulation),1);


%build the partial matrix
% for i=1:max(tri.size)
%     Triangulation(i,1)=x(tri.Triangulation);
%     Triangulation(i,2)=y(i);
%     Triangulation(i,2)=z(i);
%     
% %     centre = tri_center([Triangle_set{i}.x1 Triangle_set{i}.y1 Triangle_set{i}.z1],...
% %                                                 [Triangle_set{i}.x2 Triangle_set{i}.y2 Triangle_set{i}.z2],...
% %                                                 [Triangle_set{i}.x3 Triangle_set{i}.y3 Triangle_set{i}.z3],...
% %                                                 'circumcenter');
% %     
% end
    
[NormalVx NormalVy NormalVz PosVx PosVy PosVz]=computeNormalVectorTriangulation([x y z],tri,'center-cells');


while datenum(time, 'yyyy/mm/dd HH:MM:SS') <= datenum(tend, 'yyyy/mm/dd HH:MM:SS')
    
    [Az El] = SolarAzEl(time,lat,long,alt);
    
   
    
    %convert to rads
    A  = Az * pi/180;
    E  = El * pi/180;
    fprintf('%s\t%f\t%f\n',datestr(addtodate(datenum(time, 'yyyy/mm/dd HH:MM:SS'),-UTC_offset,'hour'), 'yyyy/mm/dd HH:MM:SS'),pi/2-E,90-El)
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
        proj_z(i) = coord(3);
    end
    
    cosi = zeros(max(size(tri)),1);


    r=1;
    sun = [cos(E)*sin(A);cos(E)*cos(A);sin(E)];

    for i=1:max(size(tri))
        angle = acos( dot(sun, [NormalVx(i)  NormalVy(i) NormalVz(i)]));
        cosi(i) = (cos(angle)+1.0)/2.0;
    end



%testing code that plots the rectangles
% ------------------------------------------------'
%     tbb=[bbx(1:4),bby(1:4)]; %ignore the last point because it is out of
%     order. Indecies follow the above naming convention
% Plots the sub rect points
%     for i=1:4
%         hold on
%        plot(tbb(i,1),tbb(i,2),'o');text(tbb(i,1)+5,tbb(i,2)+5,num2str(i)); 
%        
%     end

    %plot bounding rectangle
%     plot(bbx(:),bby(:),'color','red','linewidth',5);
% 
% % Plots the sub rects
%     hold on
%         for i=1:sg_nx
%            plot(rectangles{i}.vertex(:,1),rectangles{i}.vertex(:,2),'color','blue','linewidth',5);
%         end
%     hold off
%     pause
        
%----------------------------------------
     proj_z = (proj_z)/max(abs(proj_z)) +1 ;
    
    if  exist('p','var')==0
        if strcmp(viewpoint,'basin')
            p = patch( ...
            'Vertices',[x(:) y(:) z(:)], ...
            'Faces', tri.Triangulation, ...
            'facevertexcdata',cosi(:),...
            'facecolor','flat',...
            'edgecolor','none');
        else
    
%         for "as the sun"
            p = patch( ...
                'Vertices',[proj_x(:) proj_y(:)], ...
                'Faces', tri.Triangulation, ...
                'facevertexcdata',proj_z(:),...
                'facecolor','flat',...
                'edgecolor','none');
        end

    else
        if strcmp(viewpoint,'basin')
             set(p,'facevertexcdata',cosi(:));
        else
            %for "as the sun"
            set(p,'Vertices',[proj_x(:) proj_y(:) ],'facevertexcdata',proj_z(:));
  
            axis tight

        end
        refreshdata
    end
    if exist('ht','var')==1
        delete(ht.th);
    end

    
 
   
    if strcmp(viewpoint,'basin')
%        set(get(gca,'YLabel'),'string','Northing (m)','Fontsize',20)
%         set(get(gca,'XLabel'),'string','Easting (m)','Fontsize',20)
%         set(get(gca,'ZLabel'),'string','Elevation (m)','Fontsize',20)
%         set(gca,'Fontsize',12,'FontWeight','bold') 
%         cb=colorbar%('location','southoutside')
%         set(cb,'fontsize',16);
%         set(gcf,'color','white');
        set(gca,'visible','off');
        set(gcf,'color','black');
    else
        set(gca,'visible','off');
    end
    ht= mtit(datestr(addtodate(datenum(time, 'yyyy/mm/dd HH:MM:SS'),-UTC_offset,'hour'), 'yyyy/mm/dd HH:MM:SS'), 'fontsize',14) ;
     set(ht.th,'color','white');
    
%     F(frame)=getframe(gcf);
%      export_fig( num2str(frame),'-png');
%     frame=frame+1;
     
    
    
    time = datestr(addtodate(datenum(time,'yyyy/mm/dd HH:MM:SS'),step_size,step_type),'yyyy/mm/dd HH:MM:SS');
 
% pause

end
hold off
end


