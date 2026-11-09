function result=telemac_boundary(slfPath,cliPath,segments)
% TELEMAC_BOUNDARY opens an existing 2-D SELAFIN and exports CLI + Blue Kenue BC2.
% telemac_boundary                     Choose an SLF and open the editor.
% telemac_boundary(slfPath)            Open the editor for a known SLF.
% telemac_boundary(slfPath,cliPath,segments) exports without a GUI.
% Segment fields: Name, NodeIDs, Codes [LIHBOR LIUBOR LIVBOR LITBOR],
% Values [HBOR UBOR VBOR AUBOR TBOR ATBOR BTBOR].
if nargin==0
    [file,folder]=uigetfile('*.slf','Select an existing TELEMAC SLF file');
    if isequal(file,0), result=[]; return; end
    slfPath=fullfile(folder,file);
end
[nodes,triangles,ipobo]=readSlfGeometry(slfPath);
[ids,loops,ipobo]=boundaryLoops(nodes,triangles,ipobo);
if nargin<2
    [cliPath,bc2Path,segments]=editBoundaries(nodes,triangles,ipobo,ids,loops,slfPath);
else
    if nargin<3, segments=[]; end
    [cliPath,bc2Path]=saveBoundaries(cliPath,nodes,triangles,ipobo,ids,loops,segments,false);
end
result=struct('GeometryFile',slfPath,'BoundaryFile',cliPath,'BC2File',bc2Path,'Segments',segments);
if ~isempty(cliPath)
    fprintf('CLI: %s\nBC2: %s\n',cliPath,bc2Path);
end
end

function [xy,tri,ipobo]=readSlfGeometry(path)
fid=fopen(path,'r','ieee-be');
assert(fid>=0,'Cannot open SLF file: %s',path);
first=fread(fid,1,'int32');
if isequal(first,80)
    frewind(fid);
else
    fclose(fid);
    fid=fopen(path,'r','ieee-le');
    assert(fid>=0,'Cannot open SLF file: %s',path);
end
guard=onCleanup(@() fclose(fid));
slfRecord(fid,80,'uint8');
counts=slfRecord(fid,2,'int32');
assert(all(counts>=0),'Invalid SLF variable count.');
for i=1:sum(counts)
    slfRecord(fid,32,'uint8');
end
iparam=slfRecord(fid,10,'int32');
if iparam(10)==1
    slfRecord(fid,6,'int32');
end
dims=slfRecord(fid,4,'int32');
assert(dims(1)>0 && dims(2)>=3 && dims(3)==3 && iparam(7)<=1, ...
    'Only 2D triangular SLF geometry is supported.');
tri=reshape(slfRecord(fid,3*dims(1),'int32'),3,[]).';
ipobo=slfRecord(fid,dims(2),'int32');
bytes=fread(fid,1,'int32');
assert(~isempty(bytes) && any(bytes==dims(2)*[4 8]),'Invalid SLF coordinate record.');
assert(fseek(fid,-4,'cof')==0,'Cannot read SLF coordinate record.');
kind='single';
if bytes==8*dims(2), kind='double'; end
x=slfRecord(fid,dims(2),kind);
y=slfRecord(fid,dims(2),kind);
xy=[x y]+iparam(3:4).';
assert(all(isfinite(xy),'all') && all(tri>=1 & tri<=size(xy,1),'all'), ...
    'Invalid SLF coordinates or triangle node indices.');
end

function a=slfRecord(fid,n,kind)
bytes=4;
if strcmp(kind,'uint8'), bytes=1; end
if strcmp(kind,'double'), bytes=8; end
assert(isequal(fread(fid,1,'int32'),n*bytes),'Invalid SLF record length.');
a=fread(fid,n,kind);
tail=fread(fid,1,'int32');
assert(numel(a)==n && isequal(tail,n*bytes),'Incomplete or corrupt SLF record.');
end

function [ids,loops,rowOfNode]=boundaryLoops(xy,t,ipobo)
a=xy(t(:,2),:)-xy(t(:,1),:);
b=xy(t(:,3),:)-xy(t(:,1),:);
area2=a(:,1).*b(:,2)-a(:,2).*b(:,1);
assert(all(area2~=0),'SLF contains zero-area triangles.');
flip=area2<0;
t(flip,[2 3])=t(flip,[3 2]);
edges=[t(:,[1 2]);t(:,[2 3]);t(:,[3 1])];
[~,~,group]=unique(sort(edges,2),'rows');
count=accumarray(group,1);
assert(all(count<=2),'Invalid SLF mesh edges.');
border=edges(count(group)==1,:);
borderNodes=unique(border(:));
assert(numel(borderNodes)>=3 && size(border,1)==numel(borderNodes) && ...
    numel(unique(border(:,1)))==numel(borderNodes) && ...
    numel(unique(border(:,2)))==numel(borderNodes),'Invalid SLF boundary loops.');
nextNode=zeros(size(ipobo));
nextNode(border(:,1))=border(:,2);
ids=find(ipobo>0);
[~,order]=sort(ipobo(ids));
ids=ids(order);
if isequal(sort(ids),borderNodes) && isequal(ipobo(ids),(1:numel(ids))')
    next=ipobo(nextNode(ids));
    ends=find(next~=(2:numel(ids)+1)');
    starts=[1;ends(1:end-1)+1];
    if ~isempty(ends) && ends(end)==numel(ids) && ...
            all(next(ends)==starts) && all(ends-starts>=2)
        loops=arrayfun(@(a,b) a:b,starts,ends,'UniformOutput',false);
        rowOfNode=ipobo;
        return;
    end
end

% Some SLFs store global node IDs in IPOBO instead of the boundary row ranks.
% Derive the oriented contours from the triangles, keeping the same node IDs.
visited=false(size(ipobo));
buffer=zeros(size(borderNodes));
contours={};
areas=[];
for first=borderNodes.'
    if visited(first), continue; end
    n=0;
    node=first;
    while node>0 && ~visited(node)
        n=n+1;
        buffer(n)=node;
        visited(node)=true;
        node=nextNode(node);
    end
    assert(node==first && n>=3,'SLF boundary loop does not close.');
    ring=buffer(1:n);
    [~,start]=min(sum(xy(ring,:),2));
    ring=ring([start:end 1:start-1]);
    v=xy(ring,:)-xy(ring(1),:);
    w=v([2:end 1],:);
    areas(end+1)=sum(v(:,1).*w(:,2)-v(:,2).*w(:,1));
    contours{end+1}=ring;
end
outer=find(areas>0);
assert(isscalar(outer) && all(areas~=0),'Exactly one outer boundary loop is required.');
contours=contours([outer find(areas<0)]);
ids=vertcat(contours{:});
rowOfNode=zeros(size(ipobo));
rowOfNode(ids)=1:numel(ids);
loops=cellfun(@(ring) rowOfNode(ring).',contours,'UniformOutput',false);
fprintf('Boundary ordering derived from the SLF triangles.\n');
end

function [path,bc2Path,segments]=editBoundaries(nodes,triangles,ipobo,ids,loops,slfPath)
path='';
bc2Path='';
segments=struct('Name',{},'NodeIDs',{},'Codes',{},'Values',{});
rowOfNode=ipobo;
[boundaryTypes,presets,tracerTypes,tracerCodes]=boundaryPresets();
fig=figure('Color','w','Position',[100 100 1200 720],'ToolBar','none');
set(fig,'Name','TELEMAC - Boundary conditions','NumberTitle','off');
ax=axes('Parent',fig,'Position',[.08 .13 .89 .76]);
hold(ax,'on');
for i=1:numel(loops)
    closed=ids(loops{i}([1:end 1]));
    plot(ax,nodes(closed,1),nodes(closed,2),'.-','Color',[.36 .32 .27], ...
        'MarkerSize',9);
end
axis(ax,'equal','tight');
xlabel(ax,'X [m]'); ylabel(ax,'Y [m]');
title(ax,'SLF boundary nodes | Default 222: closed boundary');
disableDefaultInteractivity(ax);
ax.Toolbar.Visible='off';
initialX=ax.XLim;
initialY=ax.YLim;
xy=nodes(ids,:);
picked=[];
selecting=false;
dragAnchor=[];
firstMarker=gobjects(0);
status=uicontrol(fig,'Style','text','Units','normalized','Position',[.06 .02 .88 .065], ...
    'BackgroundColor','w');
addButton=button('Add segment',.035,@addSegment);
panButton=button('Pan',.19,@(~,~) updateStatus());
set(panButton,'Style','togglebutton');
homeButton=button('Fit domain',.345,@resetView);
undoButton=button('Undo last segment',.50,@undoSegment);
saveButton=button('Save CLI + BC2',.655,@saveBoundary);
skipButton=button('Skip / Close',.81,@(~,~) close(fig));
buttons=[addButton panButton homeButton undoButton saveButton skipButton];
lines=gobjects(0,1);
set(fig,'WindowButtonDownFcn',@pressMap,'WindowButtonMotionFcn',@dragMap, ...
    'WindowButtonUpFcn',@releaseMap,'WindowScrollWheelFcn',@scrollMap, ...
    'WindowKeyPressFcn',@pressKey);
updateStatus();
uiwait(fig);
if isgraphics(fig)
    cancelSelection();
    delete(buttons(isgraphics(buttons)));
    delete(status);
    set(fig,'WindowButtonDownFcn',[],'WindowButtonMotionFcn',[], ...
        'WindowButtonUpFcn',[],'WindowKeyPressFcn',[]);
    enableDefaultInteractivity(ax);
    ax.Toolbar.Visible='on';
    set(fig,'WindowScrollWheelFcn',[]);
    title(ax,'TELEMAC boundary conditions');
end
if isempty(path)
    segments=segments([]);
end

    function h=button(label,left,callback)
        h=uicontrol(fig,'Style','pushbutton','String',label,'Units','normalized', ...
            'Position',[left .935 .145 .048],'Callback',callback);
    end

    function addSegment(~,~)
        if selecting
            cancelSelection();
            return;
        end
        panButton.Value=0;
        selecting=true;
        addButton.String='Cancel selection';
        updateStatus();
    end

    function cancelSelection()
        selecting=false;
        picked=[];
        delete(firstMarker(isgraphics(firstMarker)));
        firstMarker=gobjects(0);
        addButton.String='Add segment';
        updateStatus();
    end

    function updateStatus()
        if panButton.Value
            hint='Pan: drag with the left mouse button. Toggle Pan off to resume selection.';
        elseif selecting && isempty(picked)
            hint='Click the first endpoint node.';
        elseif selecting
            hint=sprintf('First endpoint N=%d (X=%.2f, Y=%.2f). Click the second endpoint. Esc: cancel.', ...
                ids(picked),xy(picked,1),xy(picked,2));
        else
            hint='Click Add segment to define a boundary segment.';
        end
        set(status,'String',[hint '  |  Mouse wheel: zoom. Right-drag: pan.']);
    end

    function tf=pointerOnMap()
        hit=hittest(fig);
        tf=isequal(hit,ax) || isequal(ancestor(hit,'axes'),ax);
    end

    function point=mapPoint()
        point=ax.CurrentPoint(1,1:2);
    end

    function pressMap(~,~)
        if ~pointerOnMap(), return; end
        click=fig.SelectionType;
        point=mapPoint();
        if strcmp(click,'alt') || (panButton.Value && strcmp(click,'normal'))
            dragAnchor=point;
        elseif selecting && strcmp(click,'normal')
            [~,row]=min(sum((xy-point).^2,2));
            if isempty(picked)
                picked=row;
                firstMarker=plot(ax,xy(row,1),xy(row,2),'mo','MarkerFaceColor','m', ...
                    'MarkerSize',7,'HitTest','off');
                updateStatus();
            elseif picked==row
                set(status,'String','Same node selected. Click a different second endpoint node.');
            else
                ends=[picked row];
                cancelSelection();
                completeSegment(ends);
            end
        end
    end

    function dragMap(~,~)
        if isempty(dragAnchor), return; end
        delta=dragAnchor-mapPoint();
        xlim(ax,ax.XLim+delta(1));
        ylim(ax,ax.YLim+delta(2));
    end

    function releaseMap(~,~)
        dragAnchor=[];
    end

    function scrollMap(~,event)
        if ~pointerOnMap() || ~isempty(dragAnchor), return; end
        center=mapPoint();
        factor=1.2^max(-8,min(8,event.VerticalScrollCount));
        xx=center(1)+factor*(ax.XLim-center(1));
        yy=center(2)+factor*(ax.YLim-center(2));
        if all(isfinite([xx yy])) && all([diff(xx) diff(yy)]>1e-6)
            xlim(ax,xx);
            ylim(ax,yy);
        end
    end

    function resetView(~,~)
        xlim(ax,initialX);
        ylim(ax,initialY);
    end

    function pressKey(~,event)
        if strcmp(event.Key,'escape') && selecting
            cancelSelection();
        end
    end

    function completeSegment(ends)
        set(buttons,'Enable','off');
        unlock=onCleanup(@() set(buttons(isgraphics(buttons)),'Enable','on'));
        highlight=gobjects(0);
        try
            assert(ends(1)~=ends(2),'Select two different boundary nodes.');
            loop=find(cellfun(@(r) all(ismember(ends,r)),loops),1);
            assert(~isempty(loop),'Both endpoints must be on the same outer or island boundary.');
            r=loops{loop};
            a=find(r==ends(1),1); b=find(r==ends(2),1);
            forward=r(mod(a-1+(0:mod(b-a,numel(r))),numel(r))+1);
            reverse=r(mod(a-1-(0:mod(a-b,numel(r))),numel(r))+1);
            candidates={forward,reverse};
            lengths=cellfun(@(v) sum(vecnorm(diff(xy(v,:)),2,2)),candidates);
            [~,choice]=min(lengths);
            while true
                rows=candidates{choice};
                highlight=plot(ax,xy(rows,1),xy(rows,2),'m.-','LineWidth',2,'MarkerSize',12);
                drawnow;
                selection=questdlg(sprintf('Endpoints N=%d and N=%d. %d nodes, %.1f m. Select this path?', ...
                    ids(ends(1)),ids(ends(2)),numel(rows),lengths(choice)),'Confirm boundary segment','Select','Other path','Cancel','Select');
                if strcmp(selection,'Select'), break; end
                delete(highlight);
                if ~strcmp(selection,'Other path'), return; end
                choice=3-choice;
            end
            [kind,ok]=listdlg('Name','Boundary Code','PromptString','Select the boundary type', ...
                'SelectionMode','single','ListString',boundaryTypes,'ListSize',[390 170]);
            if ~ok, delete(highlight); return; end
            codes=presets(kind,:);
            if kind==numel(boundaryTypes)
                entered=inputdlg({'LIHBOR (0-6)','LIUBOR (0-6)','LIVBOR (0-6)'}, ...
                    'Custom boundary codes',1,{'2','2','2'});
                if isempty(entered), delete(highlight); return; end
                codes=cellfun(@str2double,entered).';
            end
            tracerDefault=1;
            if ~all(codes==2), tracerDefault=3; end
            [tracer,ok]=listdlg('Name','Tracer Code','PromptString','Select the tracer boundary type', ...
                'SelectionMode','single','ListString',tracerTypes,'InitialValue',tracerDefault, ...
                'ListSize',[330 90]);
            if ~ok, delete(highlight); return; end
            labels={'Segment name (use underscores instead of spaces)','HBOR','UBOR','VBOR','AUBOR','TBOR','ATBOR','BTBOR'};
            values=inputdlg(labels,'Boundary values (entire segment)',1, ...
                [{sprintf('boundary_%d',numel(segments)+1)} repmat({'0'},1,7)]);
            if isempty(values), delete(highlight); return; end
            entry=struct('Name',strtrim(values{1}),'NodeIDs',ids(rows), ...
                'Codes',[codes tracerCodes(tracer)],'Values',cellfun(@str2double,values(2:end)).');
            candidate=[segments entry];
            boundaryTable(ids,rowOfNode,loops,candidate);
            segments=candidate;
            set(highlight,'Color',[.02 .44 .78],'Marker','none','DisplayName',entry.Name);
            lines(end+1)=highlight;
            set(status,'String',sprintf('%s | %d nodes | Code: %d%d%d | Tracer: %d', ...
                entry.Name,numel(rows),entry.Codes));
        catch err
            if isgraphics(fig)
                delete(highlight(isgraphics(highlight)));
                errordlg(err.message,'Boundary selection','modal');
            end
        end
    end

    function undoSegment(~,~)
        if isempty(segments), return; end
        delete(lines(end));
        lines(end)=[];
        segments(end)=[];
        set(status,'String',sprintf('%d segments ready to save.',numel(segments)));
    end

    function saveBoundary(~,~)
        if selecting
            set(status,'String','Select the second endpoint or press Esc to cancel before saving.');
            return;
        end
        [folder,stem]=fileparts(slfPath);
        [f,p]=uiputfile('*.cli','Save boundary conditions file',fullfile(folder,[stem '.cli']));
        if isequal(f,0), return; end
        target=fullfile(p,f);
        [outFolder,outStem]=fileparts(target);
        if isfile(fullfile(outFolder,[outStem '.bc2']))
            answer=questdlg('Replace the existing BC2 file with the same name?', ...
                'Save BC2','Replace','Cancel','Cancel');
            if ~strcmp(answer,'Replace'), return; end
        end
        try
            [path,bc2Path]=saveBoundaries(target,nodes,triangles,ipobo,ids,loops,segments,true);
            uiresume(fig);
        catch err
            errordlg(err.message,'Boundary file write error','modal');
        end
    end
end

function [labels,codes,tracerLabels,tracerCodes]=boundaryPresets()
labels={'Closed boundary / wall (222)','Prescribed discharge Q (455)','Prescribed depth H (544)', ...
    'Prescribed discharge Q and depth H (555)','Prescribed velocity U and V (466)','Prescribed velocity U, V and depth H (566)', ...
    'Incident waves / ARTEMIS (111)','Custom code'};
codes=[2 2 2;4 5 5;5 4 4;5 5 5;4 6 6;5 6 6;1 1 1;0 0 0];
tracerLabels={'Closed boundary / wall (2)','Prescribed tracer (5)','Free tracer (4)'};
tracerCodes=[2 5 4];
end

function [table,names,sections]=boundaryTable(ids,rowOfNode,loops,segments)
% CONLIM columns: three hydraulic flags, four values, tracer flag,
% three tracer values, global node number and IPOBO row number.
table=zeros(numel(ids),13);
table(:,[1 2 3 8])=2;
table(:,12)=ids;
table(:,13)=(1:numel(ids))';
assigned=false(numel(ids),1);
names=strings(numel(ids),1);
sections=cell(numel(segments),1);
for i=1:numel(segments)
    s=segments(i);
    assert(~isempty(regexp(s.Name,'^[A-Za-z0-9_.-]+$','once')), ...
        'Segment names must use letters, digits, underscores, hyphens or periods.');
    validateattributes(s.NodeIDs,{'numeric'},{'vector','integer','finite','positive','<=',numel(rowOfNode)});
    validateattributes(s.Codes,{'numeric'},{'vector','numel',4,'integer','finite','>=',0,'<=',6});
    validateattributes(s.Values,{'numeric'},{'vector','numel',7,'real','finite'});
    assert(ismember(s.Codes(4),[2 4 5]),'Tracer code must be 2, 4 or 5.');
    rows=rowOfNode(s.NodeIDs(:));
    assert(numel(rows)>=2 && all(rows>0) && numel(unique(rows))==numel(rows), ...
        'A segment must contain at least two distinct boundary nodes.');
    loop=find(cellfun(@(r) all(ismember(rows,r)),loops),1);
    assert(~isempty(loop),'A segment must stay on a single boundary loop.');
    step=mod(diff(rows),numel(loops{loop}));
    assert(all(step==1) || all(step==numel(loops{loop})-1), ...
        'Segment nodes must be ordered and contiguous along the boundary.');
    assert(~any(assigned(rows)),'Segment overlaps an assigned node. Undo the last segment if needed.');
    table(rows,[1 2 3 8])=repmat(s.Codes(:).',numel(rows),1);
    table(rows,[4 5 6 7 9 10 11])=repmat(s.Values(:).',numel(rows),1);
    assigned(rows)=true;
    names(rows)=string(s.Name);
    ordered=sort(rows);
    breaks=find(diff(ordered)>1);
    ranges=[ordered([1;breaks+1]) ordered([breaks;numel(ordered)])];
    if size(ranges,1)==2
        ranges=flipud(ranges);
    end
    sections{i}=ranges;
end
end

function [cliPath,bc2Path]=saveBoundaries(cliPath,nodes,triangles,ipobo,ids,loops,segments,overwrite)
[folder,stem,ext]=fileparts(cliPath);
assert(strcmpi(ext,'.cli'),'Output file extension must be .cli.');
if isempty(folder), folder=pwd; end
cliPath=fullfile(folder,[stem ext]);
bc2Path=fullfile(folder,[stem '.bc2']);
paths={cliPath,bc2Path};
exists=cellfun(@isfile,paths);
assert(overwrite || ~any(exists),'A CLI or BC2 file with this name already exists.');
[table,names,sections]=boundaryTable(ids,ipobo,loops,segments);
stage=tempname(folder);
[ok,message]=mkdir(stage);
assert(ok,'%s',message);
guard=onCleanup(@() rmdir(stage,'s'));
staged={fullfile(stage,'out.cli'),fullfile(stage,'out.bc2')};
writeCli(staged{1},table,names);
writeBc2(staged{2},nodes,triangles,table,names,loops,sections,segments,stem);
verifyCli(staged{1},table,names);
% Prepare both files before moving either one to its final path.
backup={fullfile(stage,'previous.cli'),fullfile(stage,'previous.bc2')};
for i=find(exists)
    [ok,message]=copyfile(paths{i},backup{i});
    assert(ok,'%s',message);
end
published=false(1,2);
try
    for i=1:2
        [ok,message]=movefile(staged{i},paths{i},'f');
        assert(ok,'%s',message);
        published(i)=true;
    end
catch err
    for i=find(published)
        if exists(i)
            [ok,message]=movefile(backup{i},paths{i},'f');
            assert(ok,'%s',message);
        else
            delete(paths{i});
        end
    end
    rethrow(err);
end
end

function writeCli(path,table,names)
fid=fopen(path,'w');
assert(fid>=0,'Cannot open CLI file.');
guard=onCleanup(@() fclose(fid));
fmt='%d\t%d\t%d\t%.17g\t%.17g\t%.17g\t%.17g\t%d\t%.17g\t%.17g\t%.17g\t%d\t%d';
for i=1:size(table,1)
    assert(fprintf(fid,fmt,table(i,:))>0,'Cannot write CLI row.');
    if strlength(names(i))>0
        assert(fprintf(fid,'\t#\t%s',names(i))>0,'Cannot write CLI segment name.');
    end
    assert(fprintf(fid,'\n')>0,'Cannot write CLI row.');
end
end

function verifyCli(path,table,names)
fid=fopen(path,'r');
assert(fid>=0,'Cannot verify CLI file.');
guard=onCleanup(@() fclose(fid));
for i=1:size(table,1)
    row=fgetl(fid);
    assert(ischar(row),'Missing CLI row.');
    parts=split(string(row),'#');
    name="";
    if numel(parts)==2, name=strtrim(parts(2)); end
    assert(numel(parts)<=2 && isequal(sscanf(parts(1),'%f').',table(i,:)) && ...
        name==names(i),'CLI values or segment name do not match.');
end
assert(~ischar(fgetl(fid)),'Unexpected extra row in CLI file.');
end

function writeBc2(path,nodes,triangles,table,names,loops,sections,segments,titleText)
fid=fopen(path,'w');
assert(fid>=0,'Cannot open BC2 file.');
guard=onCleanup(@() fclose(fid));
fprintf(fid,'#########################################################################\n');
fprintf(fid,':FileType bc2  ASCII  EnSim 1.0\n# DataType Conlim Boundary Conditions Generator\n');
fprintf(fid,':Name %s\n:Title %s\n',titleText,titleText);
fprintf(fid,':NodeCount %d\n:ElementCount %d\n:ElementType T3\n',size(nodes,1),size(triangles,1));
fprintf(fid,':BoundarySegmentCount %d\n',numel(segments));
fprintf(fid,'# id code sectionCount startNode1 endNode1 startNode2 endNode2 tracerCode name\n');
for i=1:numel(segments)
    s=segments(i);
    ranges=sections{i};
    n=size(ranges,1);
    if n==1, ranges(2,:)=[1 1]; end
    fprintf(fid,':BoundarySegment %d %d%d%d %d %d %d %d %d %d "%s"\n', ...
        i,s.Codes(1),s.Codes(2),s.Codes(3),n,ranges(1,1),ranges(1,2), ...
        ranges(2,1),ranges(2,2),s.Codes(4),s.Name);
end
fprintf(fid,':ShorelineCount %d\n:ShorelineNodeCount %d\n:EndHeader\n',numel(loops),size(table,1));
fprintf(fid,':BeginNodes %d\n',size(nodes,1));
fprintf(fid,'%.17g %.17g 0\n',nodes.');
fprintf(fid,':EndNodes\n:BeginElements %d\n',size(triangles,1));
fprintf(fid,'%d %d %d\n',triangles.');
fprintf(fid,':EndElements\n:BeginTable %d 15\n',size(table,1));
shore=zeros(size(table,1),1);
for i=1:numel(loops)
    shore(loops{i})=i;
end
fmt='%d %d %d %.17g %.17g %.17g %.17g %d %.17g %.17g %.17g %d %d %d';
for i=1:size(table,1)
    fprintf(fid,fmt,[table(i,:) shore(i)]);
    fprintf(fid,'  %s\n',names(i));
end
fprintf(fid,':EndTable\n');
[message,number]=ferror(fid);
assert(number==0,'%s',message);
end
