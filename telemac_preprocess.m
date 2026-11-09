function result = telemac_preprocess(cfg)
% Saves only the requested SLF. ShowFigure=true displays one mesh/bottom figure.
% Gmsh scratch files are removed on success or error.
% Closed inner rings and an optional IslandFile define islands.
% Set IslandSpacing to 0 to disable island-coast refinement.
% Interactive mode selects rectangles, circles or polygons for refinement.
% RefinementRegions: struct array with Bounds=[xmin xmax ymin ymax],
% TargetSize (metres), Transition (metres outside the shape).
% Optional Shape: rectangle (default), circle (Center/Radius), polygon (Vertices).

interactive = nargin == 0;
if interactive, cfg = struct; end
base = fileparts(mfilename('fullpath'));
defaults = struct('OutlineFile','','BathymetryFile','', ...
    'OutputFile',fullfile(base,'geometry.slf'),'BoundarySpacing',200, ...
    'InteriorSpacing',100,'GrowthRatio',1.2,'GeometryTolerance',0.5, ...
    'GridStep',5,'SmoothIterations',10,'ZScale',1,'ZOffset',0, ...
    'MissingPolicy','error','GmshPath','','ShowFigure',true,'Overwrite',false, ...
    'IslandFile','','IslandSpacing',100,'RefinementRegions',[]);
assert(isstruct(cfg) && isscalar(cfg),'cfg must be a scalar struct.');
names = fieldnames(cfg);
for i=1:numel(names)
    assert(isfield(defaults,names{i}),'Unknown option: %s',names{i});
    defaults.(names{i}) = cfg.(names{i});
end
cfg = defaults;
if interactive, cfg = askSettings(cfg); end
validateSettings(cfg);
cfg.GmshPath = findGmsh(cfg.GmshPath);
[outdir,stem,ext] = fileparts(cfg.OutputFile);
assert(strcmpi(ext,'.slf'),'OutputFile must end in .slf');
if isempty(outdir), outdir=pwd; end
if ~isfolder(outdir), mkdir(outdir); end
cfg.OutputFile = fullfile(outdir,[stem ext]);
if ~cfg.Overwrite
    assert(~isfile(cfg.OutputFile), ...
        'Output already exists. Choose another filename or set Overwrite=true.');
end
[~,scratchName] = fileparts(tempname);
workdir = fullfile(tempdir,['telemac_' scratchName]);
mkdir(workdir);
cleanupWork = onCleanup(@() removeWorkdir(workdir));
fprintf('\n1/6 Reading outline and bathymetry...\n');
rings = readOutlines(cfg.OutlineFile);
[~,outer] = max(cellfun(@(p) abs(signedArea(p)),rings));
boundary = rings{outer};
rings(outer) = [];
if ~isempty(cfg.IslandFile)
    rings = [rings; readOutlines(cfg.IslandFile)];
end
islands = cellfun(@flipud,rings,'UniformOutput',false);
checkDomain(boundary,islands);
rawBoundary = boundary;
rawIslands = islands;
fprintf('   %d island boundaries detected.\n',numel(islands));
data = readmatrix(cfg.BathymetryFile,'FileType','text');
assert(size(data,2)==3 && size(data,1)>=3 && all(isfinite(data(:))), ...
    'Bathymetry must contain exactly three finite numeric columns: X Y Z.');
data(:,3) = cfg.ZScale*data(:,3) + cfg.ZOffset;
[xy,~,groups] = unique(data(:,1:2),'rows');
z = accumarray(groups,data(:,3),[],@mean);
duplicateCount = size(data,1)-size(xy,1);
if duplicateCount>0
    warning('Repeated bathymetry XY coordinates: %d; Z values averaged.',duplicateCount);
end
assert(rank(xy-mean(xy,1))==2,'Bathymetry points must span a 2-D area.');
assert(any(diff(boundary(:,1))) && any(diff(boundary(:,2))), ...
    'Outline must enclose a two-dimensional area.');
assert(max(max(boundary)-min(boundary))>1 && ...
    max(abs(boundary(:)))>180,'Provide projected XY coordinates in metres.');
fprintf('2/6 Resampling outer boundary and generating mesh with Gmsh...\n');
fprintf('   Requested boundary %.6g m; interior %.6g m; growth %.6g.\n', ...
    cfg.BoundarySpacing,cfg.InteriorSpacing,cfg.GrowthRatio);
if ~isempty(islands) && cfg.IslandSpacing>0
    fprintf('   Island refinement target: %.6g m.\n',min(cfg.IslandSpacing,cfg.BoundarySpacing));
end
for regionIndex=1:numel(cfg.RefinementRegions)
    r=cfg.RefinementRegions(regionIndex);
    fprintf('   Local region %d: target %.6g m; transition %.6g m.\n', ...
        regionIndex,r.TargetSize,r.Transition);
end
boundary = simplifyClosed(boundary,cfg.GeometryTolerance);
boundary = subdivideBoundary(boundary,cfg.BoundarySpacing);
for i=1:numel(islands)
    islands{i} = simplifyClosed(islands{i},cfg.GeometryTolerance);
    islands{i} = subdivideBoundary(islands{i},min(100,cfg.BoundarySpacing));
end
checkDomain(boundary,islands);
geo = fullfile(workdir,'mesh.geo');
msh = fullfile(workdir,'mesh.msh');
origin = writeGeo(geo,boundary,islands,cfg);
cmd = sprintf('%s %s -2 -format msh2 -o %s -nt 1 -v 3', ...
    shellQuote(cfg.GmshPath),shellQuote(geo),shellQuote(msh));
[status,logtext] = system(cmd);
assert(status==0 && isfile(msh),'Gmsh failed:\n%s',logtext);
[nodes,triangles] = readMsh(msh,origin);
[triangles,ipobo,quality,meshStats] = checkMesh(nodes,triangles,boundary,islands);
fprintf('   %d nodes, %d triangles, %d boundary nodes\n', ...
    size(nodes,1),size(triangles,1),nnz(ipobo));
fprintf('   Actual edge lengths: median %.6g m; maximum %.6g m.\n', ...
    meshStats.MedianEdge,meshStats.MaximumEdge);
fprintf('3/6 Interpolating bottom elevation...\n');
F = scatteredInterpolant(xy(:,1),xy(:,2),z,'natural','none');
effectiveGridStep = cfg.GridStep;
if cfg.GridStep>0
    % Keep the two-stage interpolation used in grid.m, but derive extents.
    h = cfg.GridStep;
    lo = floor(min(nodes)/h)*h-h;
    hi = ceil(max(nodes)/h)*h+h;
    cellCount = prod(floor((hi-lo)/h)+1);
    if cellCount>2e7
        fprintf('   Grid would have %.3g cells. Interpolating directly at mesh nodes.\n',cellCount);
        effectiveGridStep = 0;
    else
        gx = lo(1):h:hi(1); gy = lo(2):h:hi(2);
        [X,Y] = meshgrid(gx,gy);
        Z = F(X,Y);
        bottom = interp2(X,Y,Z,nodes(:,1),nodes(:,2),'linear',NaN);
    end
end
if effectiveGridStep==0
    bottom = F(nodes(:,1),nodes(:,2));
end
missing = ~isfinite(bottom);
missingCount = nnz(missing);
if missingCount>0
    if strcmpi(cfg.MissingPolicy,'nearest')
        Fn = scatteredInterpolant(xy(:,1),xy(:,2),z,'nearest','nearest');
        bottom(missing) = Fn(nodes(missing,1),nodes(missing,2));
        warning('%d nodes filled by explicitly selected nearest-neighbour policy.',missingCount);
    else
        error(['%d nodes lack bathymetry. No SLF written. ' ...
            'Extend the data or explicitly select MissingPolicy=''nearest''.'],missingCount);
    end
end
assert(all(isfinite(bottom)),'Nonfinite BOTTOM values cannot be written.');
fprintf('4/6 Writing SELAFIN and reading it back...\n');
temporary = fullfile(workdir,'verified.slf');
cleanupTemp = onCleanup(@() deleteIfExists(temporary));
writeSlf(temporary,nodes,triangles,ipobo,bottom);
verifySlf(temporary,nodes,triangles,ipobo,bottom);
[ok,message] = movefile(temporary,cfg.OutputFile,'f');
assert(ok,'%s',message);
clear cleanupTemp
fprintf('5/6 Preparing mesh result and on-screen preview...\n');
result = struct('Config',cfg,'Nodes',nodes,'Triangles',triangles, ...
    'BOTTOM',bottom,'IPOBO',ipobo,'Quality',quality,'MeshStatistics',meshStats, ...
    'RawBoundary',rawBoundary,'MeshedBoundary',boundary, ...
    'RawIslands',{rawIslands},'MeshedIslands',{islands},'EffectiveGridStep',effectiveGridStep, ...
    'BathymetryPointCount',size(xy,1),'DuplicateBathymetryCount',duplicateCount, ...
    'MissingBathymetryCount',missingCount,'CoordinateSystem','User supplied projected metres');
if cfg.ShowFigure
    fig=[];
    try
        fig=figure('Name','TELEMAC geometry','Color','w', ...
            'Position',[100 100 1400 650]);
        tiledlayout(fig,1,2,'TileSpacing','compact');
        nexttile; triplot(triangles,nodes(:,1),nodes(:,2),'Color',[.2 .35 .55]);
        hold on; plot(rawBoundary(:,1),rawBoundary(:,2),'r-');
        for i=1:numel(rawIslands)
            plot(rawIslands{i}(:,1),rawIslands{i}(:,2),'r-');
        end
        axis equal tight; xlabel('X [m]'); ylabel('Y [m]');
        title(sprintf('Mesh: %d nodes / %d triangles',size(nodes,1),size(triangles,1)));
        nexttile; patch('Faces',triangles,'Vertices',[nodes bottom], ...
            'FaceVertexCData',bottom,'FaceColor','interp','EdgeColor','none');
        view(2); axis equal tight; colorbar; xlabel('X [m]'); ylabel('Y [m]');
        title('BOTTOM [m]'); colormap(fig,parula);
    catch err
        if ~isempty(fig) && isgraphics(fig), close(fig); end
        warning('SLF is verified, but on-screen preview failed: %s',err.message);
    end
end
clear cleanupWork
fprintf('6/6 Complete: %s\n',cfg.OutputFile);
fprintf('   Bottom: %.3f .. %.3f m; minimum angle %.2f degrees\n', ...
    min(bottom),max(bottom),meshStats.MinimumAngle);
end

function c=askSettings(c)
candidate=fullfile(fileparts(mfilename('fullpath')),'inputs');
if ~isfolder(candidate)
    candidate=fullfile(getenv('HOME'),'Desktop','ai','matlab','inputs');
end
if ~isfile(fullfile(candidate,'outline.txt')) && isfolder(fullfile(candidate,'matlab'))
    candidate=fullfile(candidate,'matlab');
end
if ~isfolder(candidate), candidate=pwd; end
[f,p]=uigetfile({'*.txt;*.i2s;*.xyz','Outline (TXT / I2S / XYZ)'}, ...
    'Select the ordered, closed outer boundary',fullfile(candidate,'outline.txt'));
assert(~isequal(f,0),'Cancelled.'); c.OutlineFile=fullfile(p,f);
[f,p]=uigetfile({'*.txt;*.xyz;*.csv','Bathymetry (X Y Z)'}, ...
    'Select bathymetry XYZ data',fullfile(p,'data.txt'));
assert(~isequal(f,0),'Cancelled.'); c.BathymetryFile=fullfile(p,f);
choice=questdlg('Do you have an island boundary file?', ...
    'Optional island file','Select file','Skip','Skip');
if strcmp(choice,'Select file')
    [f,p]=uigetfile({'*.txt;*.i2s;*.xyz','Island boundaries (TXT / I2S / XYZ)'}, ...
        'Select the closed island boundary file',fileparts(c.OutlineFile));
    if ~isequal(f,0)
        c.IslandFile=fullfile(p,f);
    end
end
while true
a=inputdlg({'Outer boundary target point spacing [m]', ...
    'Interior target edge length [m]', ...
    'Size growth coefficient (>1; differs from Blue Kenue)', ...
    'Boundary geometry tolerance [m] (smaller = more detail)', ...
    'Bathymetry intermediate grid spacing [m] (0 = direct)', ...
    'Smoothing iterations', ...
    'Z multiplier (elevation: 1; positive depth: -1)', ...
    'Z offset [m]'},'TELEMAC SLF settings',[1 76], ...
    {'200','100','1.2','0.5','5','10','1','0'});
assert(~isempty(a),'Cancelled.'); v=cellfun(@str2double,a);
c.BoundarySpacing=v(1); c.InteriorSpacing=v(2); c.GrowthRatio=v(3);
c.GeometryTolerance=v(4); c.GridStep=v(5); c.SmoothIterations=v(6);
c.ZScale=v(7); c.ZOffset=v(8);
if any(~isfinite(v(1:2))) || any(v(1:2)<=0)
    uiwait(errordlg('Global mesh targets must be positive and finite.','Mesh targets','modal'));
    continue;
end
break;
end
validateSettings(c);
% General mesh sizes are known before selecting local refinement regions.
rings=readOutlines(c.OutlineFile);
[~,outer]=max(cellfun(@(p) abs(signedArea(p)),rings));
boundary=rings{outer}; rings(outer)=[];
if ~isempty(c.IslandFile), rings=[rings;readOutlines(c.IslandFile)]; end
islands=cellfun(@flipud,rings,'UniformOutput',false);
checkDomain(boundary,islands);
previewData=readmatrix(c.BathymetryFile,'FileType','text');
assert(size(previewData,2)==3 && size(previewData,1)>=3 && all(isfinite(previewData(:))), ...
    'Bathymetry must contain exactly three finite numeric columns: X Y Z.');
c.RefinementRegions=askRefinementRegions(boundary,islands,previewData(:,1:2),previewData(:,3),c);
for k=1:numel(c.RefinementRegions)
    r=c.RefinementRegions(k);
    % Smoothstep has a peak derivative of 1.5. Account for it in band width.
    c.RefinementRegions(k).Transition=1.5*(max(v(1:2))-r.TargetSize)/(r.GrowthRatio-1);
end

choice=questdlg(['All inputs must share a projection and use XY coordinates in metres. ' ...
    'After applying the Z multiplier/offset, values must represent bed elevation. Continue?'], ...
    'Coordinate / elevation check','Continue','Cancel','Continue');
assert(strcmp(choice,'Continue'),'Cancelled.');
choice=questdlg('What should happen if a mesh node has no bathymetry?', ...
    'Missing data','Stop','Nearest data','Stop');
assert(~isempty(choice),'Cancelled.');
if strcmp(choice,'Nearest data'), c.MissingPolicy='nearest'; end
[f,p]=uiputfile('*.slf','Save SLF geometry',c.OutputFile);
assert(~isequal(f,0),'Cancelled.'); c.OutputFile=fullfile(p,f);
% uiputfile confirms replacement of the selected SLF.
c.Overwrite=true;
end

function validateSettings(c)
assert(isfile(c.OutlineFile),'Outline file not found: %s',c.OutlineFile);
assert(isfile(c.BathymetryFile),'Bathymetry file not found: %s',c.BathymetryFile);
if ~isempty(c.IslandFile)
    assert(isfile(c.IslandFile),'Island file not found: %s',c.IslandFile);
end
for name={'BoundarySpacing','InteriorSpacing'}
    validateattributes(c.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(c.GrowthRatio,{'numeric'},{'scalar','real','finite','>',1});
validateattributes(c.GeometryTolerance,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(c.GridStep,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(c.SmoothIterations,{'numeric'},{'scalar','integer','nonnegative','<=',100});
validateattributes(c.ZScale,{'numeric'},{'scalar','real','finite','nonzero'});
validateattributes(c.ZOffset,{'numeric'},{'scalar','real','finite'});
validateattributes(c.IslandSpacing,{'numeric'},{'scalar','real','finite','nonnegative'});
if ~isempty(c.RefinementRegions)
    assert(isstruct(c.RefinementRegions) && all(isfield(c.RefinementRegions, ...
        {'Bounds','TargetSize','Transition'})), ...
        'RefinementRegions requires Bounds, TargetSize and Transition.');
    for k=1:numel(c.RefinementRegions)
        r=c.RefinementRegions(k);
        validateattributes(r.Bounds,{'numeric'},{'vector','numel',4,'real','finite'});
        assert(r.Bounds(1)<r.Bounds(2) && r.Bounds(3)<r.Bounds(4), ...
            'Refinement region bounds must be [xmin xmax ymin ymax].');
        if isfield(r,'Shape') && ~isempty(r.Shape)
            shape=validatestring(r.Shape,{'rectangle','circle','polygon'});
            if strcmp(shape,'circle')
                assert(all(isfield(r,{'Center','Radius'})),'Circle needs Center and Radius.');
                validateattributes(r.Center,{'numeric'},{'vector','numel',2,'real','finite'});
                validateattributes(r.Radius,{'numeric'},{'scalar','real','finite','positive'});
            elseif strcmp(shape,'polygon')
                assert(isfield(r,'Vertices'),'Polygon needs Vertices.');
                validateattributes(r.Vertices,{'numeric'},{'2d','ncols',2,'real','finite','nonempty'});
                assert(size(r.Vertices,1)>=3,'Polygon needs at least three vertices.');
                checkBoundary(regionOutline(r));
            end
        end
        validateattributes(r.TargetSize,{'numeric'},{'scalar','real','finite','positive'});
        validateattributes(r.Transition,{'numeric'},{'scalar','real','finite','positive'});
        assert(r.TargetSize<max(c.BoundarySpacing,c.InteriorSpacing), ...
            'Local target must be smaller than the largest background target.');
    end
end
assert(any(strcmpi(c.MissingPolicy,{'error','nearest'})),'MissingPolicy: error or nearest.');
end

function regions=askRefinementRegions(boundary,islands,xy,z,c)
regions=struct('Bounds',{},'TargetSize',{},'Transition',{},'GrowthRatio',{}, ...
    'Shape',{},'Center',{},'Radius',{},'Vertices',{});
choice=questdlg('Would you like to refine selected regions on the map?', ...
    'Local mesh refinement','Yes','No','Cancel','No');
assert(~isempty(choice) && ~strcmp(choice,'Cancel'),'Cancelled.');
if ~strcmp(choice,'Yes'), return; end
fig=figure('Name','Mesh refinement regions','NumberTitle','off', ...
    'Color','w','Position',[100 100 1200 760],'ToolBar','none','MenuBar','none');
guard=onCleanup(@() deleteRegionFigure(fig));
ax=axes('Parent',fig,'Position',[.08 .20 .86 .67]); hold(ax,'on');
idx=unique(round(linspace(1,size(xy,1),min(size(xy,1),20000))));
scatter(ax,xy(idx,1),xy(idx,2),4,z(idx),'filled','HitTest','off');
plot(ax,boundary(:,1),boundary(:,2),'k-','LineWidth',1.2);
for k=1:numel(islands)
    p=islands{k}; fill(ax,p(:,1),p(:,2),[.85 .85 .85],'EdgeColor','k');
end
axis(ax,'equal','tight'); xlabel(ax,'X [m]'); ylabel(ax,'Y [m]');
title(ax,'Choose a shape, then click Add region to define it on the map');
disableDefaultInteractivity(ax); ax.Toolbar.Visible='off';
homeX=ax.XLim; homeY=ax.YLim;
selecting=false; picked=zeros(0,2); dragAnchor=[]; preview=gobjects(0);
shapes={}; labels={}; done=false;
shape=uicontrol(fig,'Style','popupmenu','String',{'Rectangle','Circle','Polygon'}, ...
    'Units','normalized','Position',[.04 .925 .14 .05],'Callback',@(~,~) clearSelection());
add=control('Add region',.20,@startSelection);
panButton=control('Pan',.36,@(~,~) hint()); panButton.Style='togglebutton';
control('Fit domain',.52,@home); control('Remove last region',.68,@undo);
control('Done / Continue',.84,@finish);
polygonDone=uicontrol(fig,'Style','pushbutton','String','Finish polygon (Enter)', ...
    'Units','normalized','Position',[.74 .09 .23 .045],'Callback',@(~,~) completeSelection());
status=uicontrol(fig,'Style','text','Units','normalized','Position',[.04 .015 .67 .12], ...
    'BackgroundColor','w','HorizontalAlignment','left');
set(fig,'WindowButtonDownFcn',@press,'WindowButtonMotionFcn',@motion, ...
    'WindowButtonUpFcn',@(~,~) release(),'WindowScrollWheelFcn',@scroll, ...
    'WindowKeyPressFcn',@key,'CloseRequestFcn',@cancelRun);
hint(); uiwait(fig); completed=done; deleteRegionFigure(fig); assert(completed,'Cancelled.');

    function h=control(label,left,callback)
        h=uicontrol(fig,'Style','pushbutton','String',label,'Units','normalized', ...
            'Position',[left .925 .14 .05],'Callback',callback);
    end
    function hint()
        if ~selecting
            message=sprintf('%d regions defined. Choose a shape and click Add region.',numel(regions));
        elseif shape.Value==3
            message=sprintf('%d vertices selected. Click vertices; press Enter to close. Esc: cancel.',size(picked,1));
        elseif shape.Value==2 && isempty(picked)
            message='Click the circle centre.';
        elseif shape.Value==2
            message='Click the circle edge to set its radius.';
        elseif isempty(picked)
            message='Click the first rectangle corner.';
        else
            message='Click the opposite corner.';
        end
        if panButton.Value, message=[message ' Pan enabled: drag with the left mouse button.']; end
        status.String={message,'Mouse wheel: zoom | Right-drag: pan | Magenta outline: refinement region'};
        polygonDone.Enable='off';
        if selecting && shape.Value==3 && size(picked,1)>=3, polygonDone.Enable='on'; end
    end
    function clearSelection()
        selecting=false; picked=zeros(0,2);
        delete(preview(isgraphics(preview))); preview=gobjects(0);
        add.String='Add region'; hint();
    end
    function startSelection(~,~)
        if selecting, clearSelection(); return; end
        selecting=true; panButton.Value=0; add.String='Cancel selection'; hint();
    end
    function yes=onMap()
        hit=hittest(fig); yes=isequal(hit,ax) || isequal(ancestor(hit,'axes'),ax);
    end
    function p=point()
        p=ax.CurrentPoint(1,1:2);
    end
    function press(~,~)
        if ~onMap(), return; end
        p=point(); click=fig.SelectionType;
        if strcmp(click,'alt') || (panButton.Value && strcmp(click,'normal'))
            dragAnchor=p;
        elseif selecting && strcmp(click,'normal')
            if isempty(picked) || norm(p-picked(end,:))>1e-6, picked(end+1,:)=p; end
            if shape.Value~=3 && size(picked,1)==2
                completeSelection();
            else
                drawPreview(p); hint();
            end
        end
    end
    function drawPreview(p)
        delete(preview(isgraphics(preview)));
        if isempty(picked), return; end
        if shape.Value==2
            radius=norm(p-picked(1,:));theta=linspace(0,2*pi,129)';
            v=picked(1,:)+radius*[cos(theta) sin(theta)];
        elseif shape.Value==1
            v=[picked(1,1) picked(1,2);p(1) picked(1,2);p;picked(1,1) p(2);picked(1,:)];
        else
            v=[picked;p;picked(1,:)];
        end
        preview=plot(ax,v(:,1),v(:,2),'m.-','HitTest','off');
    end
    function motion(~,~)
        if ~isempty(dragAnchor)
            delta=dragAnchor-point(); ax.XLim=ax.XLim+delta(1); ax.YLim=ax.YLim+delta(2);
        elseif selecting && ~isempty(picked) && onMap(), drawPreview(point()); end
    end
    function release()
        dragAnchor=[];
    end
    function scroll(~,event)
        if ~onMap() || ~isempty(dragAnchor), return; end
        center=point(); factor=1.2^max(-8,min(8,event.VerticalScrollCount));
        xx=center(1)+factor*(ax.XLim-center(1)); yy=center(2)+factor*(ax.YLim-center(2));
        if all(isfinite([xx yy])) && min([diff(xx) diff(yy)])>1e-6
            ax.XLim=xx; ax.YLim=yy;
        end
    end
    function key(~,event)
        if strcmp(event.Key,'escape'), clearSelection();
        elseif any(strcmp(event.Key,{'return','enter'})) && selecting && shape.Value==3
            completeSelection();
        end
    end
    function home(~,~)
        ax.XLim=homeX; ax.YLim=homeY;
    end
    function completeSelection()
        if ~selecting, return; end
        r=struct('Bounds',[],'TargetSize',0,'Transition',0,'GrowthRatio',1.1, ...
            'Shape','rectangle','Center',[],'Radius',[],'Vertices',[]);
        if shape.Value==3
            if size(picked,1)<3, status.String='Select at least three vertices.'; return; end
            v=picked;
            if norm(v(end,:)-v(1,:))<1e-6, v(end,:)=[]; end
            try
                checkBoundary([v;v(1,:)]);
            catch e
                uiwait(errordlg(e.message,'Invalid polygon','modal')); return;
            end
            r.Shape='polygon';r.Vertices=v;
            r.Bounds=[min(v(:,1)) max(v(:,1)) min(v(:,2)) max(v(:,2))];
        elseif size(picked,1)<2
            return;
        elseif shape.Value==2
            r.Shape='circle';r.Center=picked(1,:);r.Radius=norm(diff(picked));
            r.Bounds=[r.Center(1)+[-r.Radius r.Radius] r.Center(2)+[-r.Radius r.Radius]];
        else
            r.Bounds=[min(picked(:,1)) max(picked(:,1)) min(picked(:,2)) max(picked(:,2))];
        end
        if any([diff(r.Bounds(1:2)) diff(r.Bounds(3:4))]<=1e-6)
            status.String='The shape must have positive area.'; return;
        end
        clearSelection(); configure(r);
    end
    function configure(r)
        domain=polyshape(boundary(:,1),boundary(:,2),'Simplify',false);
        for j=1:numel(islands)
            p=islands{j};domain=subtract(domain,polyshape(p(:,1),p(:,2),'Simplify',false));
        end
        v=regionOutline(r);selected=polyshape(v(:,1),v(:,2));
        waterArea=area(intersect(domain,selected));
        if waterArea<=0
            uiwait(errordlg('The selected region does not intersect the water domain.','Region selection','modal'));return;
        end
        while true
            a=inputdlg({'Target triangle edge length inside the region [m]', ...
                'Transition coefficient (>1; 1.05 very gentle, 1.1 gentle)'}, ...
                'Local mesh target',[1 78],{num2str(min(c.BoundarySpacing,c.InteriorSpacing)/5),'1.1'});
            if isempty(a),return;end
            values=cellfun(@str2double,a);
            if all(isfinite(values)) && values(1)>0 && values(1)<max(c.BoundarySpacing,c.InteriorSpacing) && values(2)>1,break;end
            uiwait(errordlg('Target must be positive and below the largest global mesh target; coefficient must exceed 1.','Invalid value','modal'));
        end
        width=1.5*(max(c.BoundarySpacing,c.InteriorSpacing)-values(1))/(values(2)-1);
        answer=questdlg(sprintf(['Target: %.6g m\nTransition coefficient: %.6g\nTransition band: %.6g m\n' ...
            'Approximately %.3g triangles in the selected water area (rough estimate).\nAdd this region?'], ...
            values(1),values(2),width,2.31*waterArea/values(1)^2), ...
            'Confirm region','Add','Cancel','Add');
        if ~strcmp(answer,'Add'),return;end
        r.TargetSize=values(1);r.GrowthRatio=values(2);r.Transition=width;regions(end+1)=r;
        shapes{end+1}=plot(ax,v(:,1),v(:,2),'m-','LineWidth',1.4,'HitTest','off');
        labels{end+1}=text(ax,r.Bounds(1),r.Bounds(4),sprintf('  %d: %.6g m',numel(regions),r.TargetSize), ...
            'Color','m','FontWeight','bold','HitTest','off');hint();
    end
    function undo(~,~)
        if isempty(regions),return;end
        delete(shapes{end});delete(labels{end});shapes(end)=[];labels(end)=[];regions(end)=[];hint();
    end
    function finish(~,~)
        if selecting, status.String='Finish the selection or press Esc to cancel.';return;end
        done=true;uiresume(fig);
    end
    function cancelRun(~,~)
        done=false;uiresume(fig);
    end
end

function deleteRegionFigure(fig)
if isgraphics(fig),delete(fig);end
end

function v=regionOutline(r)
shape='rectangle';if isfield(r,'Shape') && ~isempty(r.Shape),shape=r.Shape;end
switch shape
    case 'circle'
        theta=linspace(0,2*pi,257)';v=r.Center(:).'+r.Radius*[cos(theta) sin(theta)];
    case 'polygon'
        v=r.Vertices;if ~isequal(v(1,:),v(end,:)),v=[v;v(1,:)];end
    otherwise
        b=r.Bounds;v=[b(1) b(3);b(2) b(3);b(2) b(4);b(1) b(4);b(1) b(3)];
end
end

function rings=readOutlines(path)
[~,~,ext]=fileparts(path);
rings = {};
if strcmpi(ext,'.i2s')
    lines=splitlines(string(fileread(path)));
    k=find(startsWith(strtrim(lines),':EndHeader'),1);
    assert(~isempty(k),'Invalid I2S header.');
    lines=lines(k+1:end); lines=lines(strlength(strtrim(lines))>0);
    i=1;
    while i<=numel(lines)
        header=sscanf(lines(i),'%f');
        assert(~isempty(header),'Invalid I2S polyline header.');
        n=header(1);
        assert(n>=4 && n==fix(n) && i+n<=numel(lines),'Invalid I2S polyline length.');
        p=zeros(n,2);
        for j=1:n
            v=sscanf(lines(i+j),'%f');
            assert(numel(v)==2,'I2S must contain XY pairs.');
            p(j,:)=v.';
        end
        rings{end+1,1}=normalizeRing(p);
        i=i+n+1;
    end
else
    % Read explicit NaN separators without discarding the boundary between rings.
    lines=strip(splitlines(string(fileread(path))));
    lines=lines(strlength(lines)>0);
    xy=nan(numel(lines),2);
    for i=1:numel(lines)
        fields=regexp(char(lines(i)),'[,;\s]+','split');
        assert(numel(fields)>=2 && numel(fields)<=3,'Outline must contain XY or XYZ columns.');
        if all(strcmpi(fields(1:2),'nan')), continue; end
        xy(i,:)=str2double(fields(1:2));
        assert(all(isfinite(xy(i,:))),'Invalid outline coordinates at line %d.',i);
    end
    first=1;
    while first<=size(xy,1)
        if isnan(xy(first,1))
            first=first+1;
            continue;
        end
        last=first+3;
        while last<=size(xy,1) && all(isfinite(xy(last,:))) && ...
                norm(xy(last,:)-xy(first,:))>=1e-5
            last=last+1;
        end
        assert(last<=size(xy,1) && all(isfinite(xy(first:last,:)),'all'), ...
            'Each outline ring must be ordered and explicitly closed (line %d).',first);
        rings{end+1,1}=normalizeRing(xy(first:last,:));
        first=last+1;
    end
end
assert(~isempty(rings),'No closed outlines in %s.',path);
end

function xy=normalizeRing(xy)
assert(size(xy,1)>=4 && all(isfinite(xy(:))),'Invalid outline coordinates.');
assert(norm(xy(1,:)-xy(end,:))<1e-5,'Each outline ring must be explicitly closed.');
xy(end,:)=xy(1,:);
xy=xy([true;vecnorm(diff(xy),2,2)>1e-8],:);
assert(size(unique(xy(1:end-1,:),'rows'),1)==size(xy,1)-1, ...
    'Outline repeats nonterminal vertices.');
if signedArea(xy)<0, xy=flipud(xy); end
end

function checkDomain(boundary,islands)
checkBoundary(boundary);
for i=1:numel(islands)
    p=islands{i};
    checkBoundary(p);
    [inside,on]=inpolygon(p(:,1),p(:,2),boundary(:,1),boundary(:,2));
    assert(all(inside & ~on) && ~ringsTouch(p,boundary), ...
        'Island %d must lie strictly inside the outer boundary.',i);
    for j=1:i-1
        q=islands{j};
        assert(~inpolygon(p(1,1),p(1,2),q(:,1),q(:,2)) && ...
            ~inpolygon(q(1,1),q(1,2),p(:,1),p(:,2)) && ~ringsTouch(p,q), ...
            'Islands %d and %d overlap, touch, or are nested.',j,i);
    end
end
end

function touched=ringsTouch(p,q)
touched=false;
if any(max(p)<min(q)) || any(max(q)<min(p)), return; end
c=q(1:end-1,:); d=q(2:end,:);
lo=min(c,d); hi=max(c,d);
for i=1:size(p,1)-1
    a=p(i,:); b=p(i+1,:); ab=b-a;
    use=all(max(a,b)>=lo-1e-8 & min(a,b)<=hi+1e-8,2);
    if ~any(use), continue; end
    x=c(use,:); y=d(use,:); cd=y-x;
    v1=cross2(ab,x-a); v2=cross2(ab,y-a);
    v3=cross2(cd,a-x); v4=cross2(cd,b-x);
    if any(v1.*v2<=1e-14 & v3.*v4<=1e-14)
        touched=true;
        return;
    end
end
end

function checkBoundary(p)
assert(abs(signedArea(p))>1,'Outline has negligible area.');
n=size(p,1)-1;
for i=1:n
    a=p(i,:); b=p(i+1,:);
    for j=i+2:n
        if i==1 && j==n, continue; end
        c=p(j,:); d=p(j+1,:);
        if any(max(a,b)<min(c,d)-1e-8) || any(max(c,d)<min(a,b)-1e-8), continue; end
        ab=b-a; cd=d-c;
        v=[cross2(ab,c-a),cross2(ab,d-a),cross2(cd,a-c),cross2(cd,b-c)];
        tol=1e-9*max([norm(ab),norm(cd),1]);
        assert(~(v(1)*v(2)<=tol^2 && v(3)*v(4)<=tol^2), ...
            'Outline intersects/touches itself near segments %d and %d.',i,j);
    end
end
end

function p=simplifyClosed(p,tolerance)
if tolerance==0, return; end
[~,k]=max(sum((p-p(1,:)).^2,2));
p=[rdp(p(1:k,:),tolerance);rdpTail(p(k:end,:),tolerance)];
end
function p=rdpTail(p,tolerance)
p=rdp(p,tolerance); p=p(2:end,:);
end
function p=rdp(p,tol)
if size(p,1)<=2, return; end
v=p(end,:)-p(1,:);
t=max(0,min(1,((p-p(1,:))*v.')/sum(v.^2)));
dist=vecnorm(p-(p(1,:)+t*v),2,2);
[d,k]=max(dist);
if d>tol
    a=rdp(p(1:k,:),tol); b=rdp(p(k:end,:),tol); p=[a;b(2:end,:)];
else
    p=p([1 end],:);
end
end
function q=subdivideBoundary(p,h)
q=cell(size(p,1)-1,1);
for i=1:numel(q)
    count=max(1,ceil(norm(p(i+1,:)-p(i,:))/h));
    q{i}=p(i,:)+(0:count-1)'/count*(p(i+1,:)-p(i,:));
end
q=vertcat(q{:}); q=[q;q(1,:)];
end

function path=findGmsh(path)
if ~isempty(path), assert(isfile(path),'Gmsh executable not found.'); return; end
if ismac
    candidates={'/Applications/Gmsh.app/Contents/MacOS/gmsh', ...
        '/opt/homebrew/bin/gmsh','/usr/local/bin/gmsh'};
elseif ispc
    candidates={fullfile(getenv('ProgramFiles'),'gmsh','gmsh.exe')};
else
    candidates={'/usr/bin/gmsh','/usr/local/bin/gmsh'};
end
for i=1:numel(candidates)
    if isfile(candidates{i}), path=candidates{i}; return; end
end
lookup='command -v gmsh'; if ispc, lookup='where gmsh'; end
[status,out]=system(lookup);
if status==0
    lines=splitlines(strtrim(out)); path=char(lines(1));
end
assert(~isempty(path) && isfile(path), ...
    'Install Gmsh or provide cfg.GmshPath (https://gmsh.info).');
end

function origin=writeGeo(path,p,islands,c)
fid=fopen(path,'w'); assert(fid>=0,'Cannot write GEO.');
guard=onCleanup(@() fclose(fid));
% Local coordinates improve geometric calculations. Restore absolute XY later.
origin=floor(min(p)); p=p-origin; n=size(p,1)-1;
lengths=vecnorm(diff(p),2,2);
% Short retained coastal segments and narrow channels need local refinement.
% Keep point sizes on coastal details; the background fields control the
% transition into the interior independently of boundary interpolation.
localSize=min(c.BoundarySpacing,min(lengths,circshift(lengths,1)));
fprintf(fid,'// XY origin: %.17g %.17g\n',origin);
for i=1:n
    fprintf(fid,'Point(%d) = {%.17g, %.17g, 0, %.17g};\n',i,p(i,1),p(i,2),localSize(i));
end
for i=1:n, fprintf(fid,'Line(%d) = {%d,%d};\n',i,i,mod(i,n)+1); end
fprintf(fid,'Curve Loop(1) = {1:%d};\n',n);
lastCurve=n;
islandTarget=c.BoundarySpacing;
if c.IslandSpacing>0
    islandTarget=min(islandTarget,c.IslandSpacing);
end
for k=1:numel(islands)
    ring=islands{k}-origin;
    m=size(ring,1)-1;
    lengths=vecnorm(diff(ring),2,2);
    sizes=min(islandTarget,min(lengths,circshift(lengths,1)));
    for j=1:m
        fprintf(fid,'Point(%d) = {%.17g, %.17g, 0, %.17g};\n', ...
            lastCurve+j,ring(j,1),ring(j,2),sizes(j));
    end
    for j=1:m
        fprintf(fid,'Line(%d) = {%d,%d};\n',lastCurve+j,lastCurve+j,lastCurve+mod(j,m)+1);
    end
    fprintf(fid,'Curve Loop(%d) = {%d:%d};\n',k+1,lastCurve+1,lastCurve+m);
    lastCurve=lastCurve+m;
end
fprintf(fid,'Plane Surface(1) = {1:%d};\n',1+numel(islands));
fprintf(fid,'Physical Curve(1) = {1:%d};\nPhysical Surface(2) = {1};\n',lastCurve);
fprintf(fid,'Field[1] = Distance;\nField[1].CurvesList = {1:%d};\nField[1].Sampling = 30;\n',n);
fprintf(fid,'Field[2] = MathEval;\n');
if c.InteriorSpacing>=c.BoundarySpacing
    fprintf(fid,'Field[2].F = "Min(%.17g, %.17g + %.17g*F1)";\n', ...
        c.InteriorSpacing,c.BoundarySpacing,c.GrowthRatio-1);
else
    fprintf(fid,'Field[2].F = "Max(%.17g, %.17g - %.17g*F1)";\n', ...
        c.InteriorSpacing,c.BoundarySpacing,c.GrowthRatio-1);
end
sizeFields=2;
if ~isempty(islands) && c.IslandSpacing>0
    fprintf(fid,'Field[5] = Distance;\nField[5].CurvesList = {%d:%d};\nField[5].Sampling = 30;\n',n+1,lastCurve);
    fprintf(fid,'Field[6] = MathEval;\nField[6].F = "%.17g + %.17g*F5";\n',islandTarget,c.GrowthRatio-1);
    sizeFields(end+1)=6;
end
% Independent local boxes refine water both around and between islands.
% Shape-based fields use cubic smoothstep with zero slope at both ends,
% avoiding an abrupt start or end to the size ramp.
for k=1:numel(c.RefinementRegions)
    r=c.RefinementRegions(k);
    id=10+k;
    distance=regionDistanceExpression(r,origin);
    blend=sprintf('Min(1,(%s)/%.17g)',distance,r.Transition);
    % Inline the blend: repeated nested MathEval field references can stall
    % older Gmsh releases during curve integration.
    expression=sprintf('%.17g+%.17g*(3*(%s)^2-2*(%s)^3)', ...
        r.TargetSize,max(c.InteriorSpacing,c.BoundarySpacing)-r.TargetSize,blend,blend);
    writeMathField(fid,id,expression);
    sizeFields(end+1)=id;
end
fprintf(fid,'Field[4] = Min;\nField[4].FieldsList = {%s};\nBackground Field = 4;\n', ...
    strjoin(string(sizeFields),','));
% Distance fields already define the transition from all boundaries.
% Extending boundary sizes here would cap a coarser interior target.
extendBoundary=0;
fprintf(fid,['Mesh.MeshSizeFromPoints = 1;\n' ...
    'Mesh.MeshSizeFromCurvature = 0;\nMesh.MeshSizeExtendFromBoundary = %d;\n' ...
    'Mesh.Algorithm = 6;\nMesh.Smoothing = %d;\nMesh.ElementOrder = 1;\n' ...
    'Mesh.RecombineAll = 0;\nMesh.MshFileVersion = 2.2;\nMesh.Binary = 0;\n'],extendBoundary,c.SmoothIterations);
end

function writeMathField(fid,id,expression)
fprintf(fid,'Field[%d] = MathEval;\nField[%d].F = ',id,id);
% Gmsh's GEO lexer limits individual string tokens. Join short tokens for
% polygons with many vertices instead of emitting one oversized literal.
if numel(expression)<=500
    fprintf(fid,'"%s";\n',expression);
else
    fprintf(fid,'StrCat(');
    for first=1:500:numel(expression)
        if first>1,fprintf(fid,',');end
        fprintf(fid,'"%s"',expression(first:min(first+499,end)));
    end
    fprintf(fid,');\n');
end
end

function expression=regionDistanceExpression(r,origin)
shape='rectangle';
if isfield(r,'Shape') && ~isempty(r.Shape),shape=lower(char(r.Shape));end
switch shape
    case 'circle'
        center=r.Center(:).'-origin;
        expression=sprintf('Max(0,Sqrt((x-(%.17g))^2+(y-(%.17g))^2)-%.17g)', ...
            center(1),center(2),r.Radius);
    case 'polygon'
        v=r.Vertices-origin;
        poly=polyshape(v(:,1),v(:,2),'Simplify',false);
        tr=triangulation(poly);terms=cell(size(tr.ConnectivityList,1),1);
        % Union of triangle half-plane distances: zero exactly in the polygon,
        % positive and continuous outside, including concave selections.
        for k=1:numel(terms)
            q=tr.Points(tr.ConnectivityList(k,:),:);q=[q;q(1,:)];
            if signedArea(q)<0,q=flipud(q);end
            edges=diff(q);normals=[edges(:,2) -edges(:,1)]./vecnorm(edges,2,2);
            offsets=-sum(normals.*q(1:3,:),2);side=cell(3,1);
            for j=1:3
                side{j}=sprintf('(%.17g*x+(%.17g)*y+(%.17g))',normals(j,1),normals(j,2),offsets(j));
            end
            terms{k}=sprintf('Max(0,Max(%s,Max(%s,%s)))',side{:});
        end
        expression=terms{1};
        for k=2:numel(terms),expression=sprintf('Min(%s,%s)',expression,terms{k});end
    otherwise
        b=r.Bounds(:).'-[origin(1) origin(1) origin(2) origin(2)];
        expression=sprintf('Sqrt(Max(0,Max(%.17g-x,x-%.17g))^2+Max(0,Max(%.17g-y,y-%.17g))^2)',b);
end
end

function [p,t]=readMsh(path,origin)
fid=fopen(path,'r'); assert(fid>=0,'Cannot read Gmsh mesh.');
guard=onCleanup(@() fclose(fid)); p=[]; t=[]; ids=[];
while ~feof(fid)
    line=strtrim(fgetl(fid));
    if strcmp(line,'$Nodes')
        n=str2double(fgetl(fid)); a=fscanf(fid,'%f',[4 n]).';
        ids=a(:,1); p=a(:,2:3);
    elseif strcmp(line,'$Elements')
        n=str2double(fgetl(fid)); t=zeros(n,3); k=0;
        for i=1:n
            a=sscanf(fgetl(fid),'%f').';
            if a(2)==2
                k=k+1; t(k,:)=a(4+a(3):6+a(3));
            end
        end
        t=t(1:k,:);
    end
end
assert(~isempty(p) && ~isempty(t),'Gmsh did not produce a triangular mesh.');
[found,t]=ismember(t,ids); assert(all(found(:)),'Invalid Gmsh node IDs.');
used=unique(t(:)); [~,t]=ismember(t,used); p=p(used,:);
p=p+origin;
end

function [t,ipobo,q,stats]=checkMesh(p,t,boundary,islands)
a=p(t(:,2),:)-p(t(:,1),:); b=p(t(:,3),:)-p(t(:,1),:);
area2=a(:,1).*b(:,2)-a(:,2).*b(:,1);
assert(all(abs(area2)>1e-8),'Degenerate triangles.');
flip=area2<0; t(flip,[2 3])=t(flip,[3 2]); area2=abs(area2);
assert(size(unique(sort(t,2),'rows'),1)==size(t,1),'Duplicate triangles.');
edges=[t(:,[1 2]);t(:,[2 3]);t(:,[3 1])];
[u,~,j]=unique(sort(edges,2),'rows'); count=accumarray(j,1);
assert(all(count<=2),'Nonmanifold mesh edges.');
G=graph(u(:,1),u(:,2),[],size(p,1)); assert(max(conncomp(G))==1,'Disconnected mesh.');
be=edges(count(j)==1,:);
assert(size(be,1)>=3,'Missing boundary.');
assert(all(accumarray(be(:,1),1,[size(p,1),1])<=1) && ...
    all(accumarray(be(:,2),1,[size(p,1),1])<=1),'Branching boundary.');
next=zeros(size(p,1),1); next(be(:,1))=be(:,2);
ipobo=zeros(size(p,1),1);
number=0; loopCount=0;
for start=find(next).'
    if ipobo(start)>0, continue; end
    loopCount=loopCount+1;
    node=start;
    while node>0 && ipobo(node)==0
        number=number+1;
        ipobo(node)=number;
        node=next(node);
    end
    assert(node==start,'Broken boundary loop.');
end
assert(loopCount==1+numel(islands),'Unexpected number of mesh boundary loops.');
centres=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
assert(all(inpolygon(centres(:,1),centres(:,2),boundary(:,1),boundary(:,2))), ...
    'Triangles outside the domain.');
expectedArea=abs(signedArea(boundary));
for i=1:numel(islands)
    ring=islands{i};
    near=all(centres>=min(ring) & centres<=max(ring),2);
    assert(~any(inpolygon(centres(near,1),centres(near,2),ring(:,1),ring(:,2))), ...
        'Triangles inside island %d.',i);
    expectedArea=expectedArea-abs(signedArea(ring));
end
area=sum(area2)/2;
relative=abs(area-expectedArea)/area;
assert(relative<1e-7,'Mesh area differs from the water domain.');
l1=vecnorm(p(t(:,2),:)-p(t(:,3),:),2,2);
l2=vecnorm(p(t(:,3),:)-p(t(:,1),:),2,2);
l3=vecnorm(p(t(:,1),:)-p(t(:,2),:),2,2);
q=2*sqrt(3)*area2./(l1.^2+l2.^2+l3.^2);
cosines=[(l2.^2+l3.^2-l1.^2)./(2*l2.*l3), ...
    (l1.^2+l3.^2-l2.^2)./(2*l1.*l3),(l1.^2+l2.^2-l3.^2)./(2*l1.*l2)];
minimum=min(acosd(max(-1,min(1,cosines))),[],'all');
edgeLengths=vecnorm(p(u(:,1),:)-p(u(:,2),:),2,2);
stats=struct('MinimumAngle',minimum,'MedianEdge',median(edgeLengths), ...
    'MaximumEdge',max(edgeLengths),'Area',area,'AreaRelativeError',relative);
if minimum<10, warning('Minimum triangle angle is %.2f degrees. Inspect the mesh.',minimum); end
end

function writeSlf(path,p,t,ib,z)
encoded=double(single(p));
assert(size(unique(encoded,'rows'),1)==size(p,1), ...
    'Float32 XY would merge nodes. Reduce coordinates with an origin or use SERAFIND.');
a=encoded(t(:,2),:)-encoded(t(:,1),:); b=encoded(t(:,3),:)-encoded(t(:,1),:);
assert(all(a(:,1).*b(:,2)-a(:,2).*b(:,1)>0), ...
    'Float32 XY would invert/collapse triangles. Use an origin or SERAFIND.');
fid=fopen(path,'w','ieee-be'); assert(fid>=0,'Cannot create SLF.');
guard=onCleanup(@() fclose(fid));
record(fid,uint8([pad('MATLAB TELEMAC GEOMETRY',72,'right') 'SERAFIN ']),'uint8');
record(fid,[1 0],'int32');
record(fid,uint8([pad('BOTTOM',16,'right') pad('M',16,'right')]),'uint8');
record(fid,[1 zeros(1,9)],'int32');
record(fid,[size(t,1) size(p,1) 3 1],'int32');
record(fid,reshape(t.',[],1),'int32');
record(fid,ib,'int32');
record(fid,p(:,1),'single'); record(fid,p(:,2),'single');
record(fid,0,'single'); record(fid,z,'single');
end
function record(fid,a,type)
bytes=4*numel(a); if strcmp(type,'uint8'), bytes=numel(a); end
assert(fwrite(fid,bytes,'int32')==1,'SLF write failed.');
assert(fwrite(fid,a,type)==numel(a),'SLF write failed.');
assert(fwrite(fid,bytes,'int32')==1,'SLF write failed.');
end
function verifySlf(path,p,t,ib,z)
fid=fopen(path,'r','ieee-be'); assert(fid>=0,'Cannot verify SLF.');
guard=onCleanup(@() fclose(fid));
title=readRecord(fid,'uint8'); assert(numel(title)==80,'Invalid SLF title.');
assert(isequal(readRecord(fid,'int32'),[1;0]),'Invalid variable counts.');
name=char(readRecord(fid,'uint8').'); assert(strcmp(strtrim(name(1:16)),'BOTTOM'));
assert(isequal(readRecord(fid,'int32'),[1;zeros(9,1)]),'Invalid IPARAM.');
assert(isequal(readRecord(fid,'int32'),[size(t,1);size(p,1);3;1]),'Invalid mesh counts.');
assert(isequal(readRecord(fid,'int32'),reshape(t.',[],1)),'Connectivity mismatch.');
assert(isequal(readRecord(fid,'int32'),ib),'IPOBO mismatch.');
assert(isequal(readRecord(fid,'single'),double(single(p(:,1)))),'X mismatch.');
assert(isequal(readRecord(fid,'single'),double(single(p(:,2)))),'Y mismatch.');
assert(isequal(readRecord(fid,'single'),0),'Invalid time.');
assert(isequal(readRecord(fid,'single'),double(single(z))),'BOTTOM mismatch.');
assert(isempty(fread(fid,1,'uint8')),'Unexpected data after BOTTOM.');
end
function a=readRecord(fid,type)
n=fread(fid,1,'int32'); assert(~isempty(n) && n>=0,'Missing SLF record.');
bytes=4; if strcmp(type,'uint8'), bytes=1; end
assert(mod(n,bytes)==0,'Invalid record length.');
a=fread(fid,n/bytes,['*' type]); a=double(a);
tail=fread(fid,1,'int32'); assert(isequal(n,tail),'SLF record marker mismatch.');
assert(numel(a)==n/bytes,'Truncated SLF record.');
end
function a=signedArea(p)
v=p-p(1,:); a=sum(v(1:end-1,1).*v(2:end,2)-v(2:end,1).*v(1:end-1,2))/2;
end
function a=cross2(u,v)
a=u(:,1).*v(:,2)-u(:,2).*v(:,1);
end
function removeWorkdir(path)
if isfolder(path)
    [ok,message]=rmdir(path,'s');
    if ~ok, warning('Could not remove temporary mesh files: %s',message); end
end
end
function deleteIfExists(path)
if isfile(path), delete(path); end
end
function q=shellQuote(s)
assert(~contains(s,newline) && ~contains(s,char(13)),'Invalid path.');
if ispc
    assert(~any(ismember(s,'"%%!&|<>^')),'Unsupported characters in Windows path.');
    q=['"' s '"'];
else
    q=['''' strrep(s,'''','''"''"''') ''''];
end
end
