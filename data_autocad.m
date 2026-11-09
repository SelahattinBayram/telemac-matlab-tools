function outputFile=data_autocad(inputFile,outputFile,openDrawing,layerStep)
% DATA_AUTOCAD reads XYZ bathymetry and opens it in AutoCAD as 3D points.
% Run without arguments from the MATLAB editor.
% Optional: data_autocad(inputFile,outputFile,false) exports without opening.
% Running without an input file opens an XYZ data file picker.
% POINT entities always retain original X,Y,Z.
% Up to 240 distinct Z values: one layer per Z. Larger sets: automatic Z bands.
% Optional fourth argument sets the layer band width in Z units (e.g. 10).

base=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(inputFile)
    startFolder=fullfile(base,'inputs');
    if ~isfolder(startFolder), startFolder=base; end
    [file,folder]=uigetfile({'*.txt;*.xyz;*.csv','XYZ data (TXT / XYZ / CSV)'}, ...
        'Select XYZ data for AutoCAD',startFolder);
    if isequal(file,0), outputFile=''; return; end
    inputFile=fullfile(folder,file);
end
assert(isfile(inputFile),'Input file not found: %s',inputFile);
if nargin<2 || isempty(outputFile)
    [~,name]=fileparts(inputFile);
    outputFile=fullfile(base,[name '.dxf']);
    if isfile(outputFile)
        [file,folder]=uiputfile('*.dxf','Save AutoCAD DXF file',outputFile);
        if isequal(file,0), outputFile=''; return; end
        outputFile=fullfile(folder,file);
    end
end
if nargin<3, openDrawing=true; end
validateattributes(openDrawing,{'logical','numeric'},{'scalar','binary'});
[folder,~,extension]=fileparts(outputFile);
assert(strcmpi(extension,'.dxf'),'Output extension must be .dxf.');
if isempty(folder), folder=pwd; outputFile=fullfile(folder,outputFile); end
assert(isfolder(folder),'Output folder not found: %s',folder);

xyz=readmatrix(inputFile,'FileType','text');
assert(~isempty(xyz) && size(xyz,2)==3 && all(isfinite(xyz(:))), ...
    'Input must contain exactly three finite numeric columns: X Y Z.');
if nargin<4, layerStep=[]; end
if ~isempty(layerStep)
    validateattributes(layerStep,{'numeric'},{'scalar','real','finite','positive'});
end
lo=min(xyz,[],1);
hi=max(xyz,[],1);
center=(lo+hi)/2;
aspect=1.5;
height=1.15*max([hi(2)-lo(2),(hi(1)-lo(1))/aspect,1]);

temporary=[tempname(folder) '.dxf'];
cleanup=onCleanup(@() deleteTemporary(temporary));
[layerCount,usedStep]=writePoints(temporary,xyz,lo,hi,center,height,aspect,layerStep);
[ok,message]=movefile(temporary,outputFile,'f');
assert(ok,'%s',message);
fprintf('%d 3D points, %d Z layers\nDXF: %s\n', ...
    size(xyz,1),layerCount,outputFile);
if usedStep>0
    fprintf('Layer elevation interval: %.6g (original point elevations preserved)\n',usedStep);
end
fprintf('Z range: %.6g to %.6g\n',lo(3),hi(3));
if ~openDrawing, return; end

if ismac
    apps=dir('/Applications/Autodesk/AutoCAD */AutoCAD *.app');
    app='';
    for i=1:numel(apps)
        if ~contains(apps(i).name,'Plot Style') && ~startsWith(apps(i).name,'Remove')
            app=fullfile(apps(i).folder,apps(i).name);
        end
    end
    assert(~isempty(app), ...
        'DXF created. AutoCAD not found; use File > Open in AutoCAD.');
    command=sprintf('open -a %s %s',quotePath(app),quotePath(outputFile));
    [status,message]=system(command);
    assert(status==0,'DXF created, but AutoCAD could not be opened: %s',message);
elseif ispc
    winopen(outputFile);
else
    error('DXF created. Use File > Open in AutoCAD: %s',outputFile);
end
end

function [layerCount,usedStep]=writePoints(path,xyz,lo,hi,center,height,aspect,layerStep)
% ASCII DXF R12 supports 3-D POINT entities without an additional library.
[zLevels,~,groups]=unique(xyz(:,3));
usedStep=0;
if isempty(layerStep) && numel(zLevels)>240
    % Choose a readable interval with about 200 bands at most.
    target=(hi(3)-lo(3))/200;
    decade=10^floor(log10(target));
    choices=[1 2 5 10]*decade;
    layerStep=choices(find(choices>=target,1));
end
if ~isempty(layerStep)
    usedStep=layerStep;
    bandIndex=floor(xyz(:,3)/layerStep);
    [bands,~,groups]=unique(bandIndex);
    zLevels=bands*layerStep;
    layerNames=arrayfun(@(z) sprintf('Z_%.12g_to_%.12g',z,z+layerStep), ...
        zLevels,'UniformOutput',false);
else
    layerNames=arrayfun(@(z) sprintf('Z_%.17g',z),zLevels,'UniformOutput',false);
end
layerCount=numel(zLevels);
assert(layerCount<=32766, ...
    'Too many layers. Use a larger elevation interval as the fourth argument.');
colors=10+floor((0:layerCount-1)*240/layerCount);
fid=fopen(path,'w');
assert(fid>=0,'Cannot write DXF file.');
guard=onCleanup(@() fclose(fid));
fprintf(fid,'0\nSECTION\n2\nHEADER\n9\n$ACADVER\n1\nAC1009\n');
fprintf(fid,'9\n$EXTMIN\n10\n%.17g\n20\n%.17g\n30\n%.17g\n',lo);
fprintf(fid,'9\n$EXTMAX\n10\n%.17g\n20\n%.17g\n30\n%.17g\n',hi);
fprintf(fid,'9\n$PDMODE\n70\n0\n9\n$PDSIZE\n40\n-0.5\n');
fprintf(fid,'9\n$LUNITS\n70\n2\n9\n$LUPREC\n70\n6\n0\nENDSEC\n');

fprintf(fid,'0\nSECTION\n2\nTABLES\n');
fprintf(fid,'0\nTABLE\n2\nLTYPE\n70\n1\n');
fprintf(fid,'0\nLTYPE\n2\nCONTINUOUS\n70\n0\n3\nSolid line\n72\n65\n73\n0\n40\n0\n0\nENDTAB\n');
fprintf(fid,'0\nTABLE\n2\nLAYER\n70\n%d\n',numel(zLevels)+1);
fprintf(fid,'0\nLAYER\n2\n0\n70\n0\n62\n7\n6\nCONTINUOUS\n');
for i=1:numel(zLevels)
    fprintf(fid,'0\nLAYER\n2\n%s\n70\n0\n62\n%d\n6\nCONTINUOUS\n',layerNames{i},colors(i));
end
fprintf(fid,'0\nENDTAB\n');
% Top view centered on the bathymetry avoids opening near the coordinate origin.
fprintf(fid,'0\nTABLE\n2\nVPORT\n70\n1\n0\nVPORT\n2\n*ACTIVE\n70\n0\n');
fprintf(fid,'10\n0\n20\n0\n11\n1\n21\n1\n');
fprintf(fid,'12\n%.17g\n22\n%.17g\n',center(1:2));
fprintf(fid,'13\n0\n23\n0\n14\n10\n24\n10\n15\n10\n25\n10\n');
fprintf(fid,'16\n0\n26\n0\n36\n1\n17\n0\n27\n0\n37\n0\n');
fprintf(fid,'40\n%.17g\n41\n%.17g\n42\n50\n43\n0\n44\n0\n50\n0\n51\n0\n',height,aspect);
fprintf(fid,'71\n0\n72\n100\n73\n1\n74\n3\n75\n0\n76\n0\n77\n0\n78\n0\n0\nENDTAB\n');
fprintf(fid,'0\nENDSEC\n0\nSECTION\n2\nBLOCKS\n0\nENDSEC\n');
fprintf(fid,'0\nSECTION\n2\nENTITIES\n');
for i=1:size(xyz,1)
    fprintf(fid,'0\nPOINT\n8\n%s\n62\n256\n10\n%.17g\n20\n%.17g\n30\n%.17g\n', ...
        layerNames{groups(i)},xyz(i,:));
end
fprintf(fid,'0\nENDSEC\n0\nEOF\n');
[message,number]=ferror(fid);
assert(number==0,'%s',message);
end

function deleteTemporary(path)
if isfile(path), delete(path); end
end

function q=quotePath(path)
assert(~contains(path,newline) && ~contains(path,char(13)),'Invalid file path.');
q=['''' strrep(path,'''','''"''"''') ''''];
end
