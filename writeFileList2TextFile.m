function writeFileList2TextFile(dDir,outFn,varargin)
% ** function writeFileList2TextFile(dDir,outFn,varargin)
% generates a list of files in a directory including subdirectories and
% saves the full file names including pathsinto a text file, one entry per
% line. Files can be selected according to their extension, and can be
% sorted according to dates (dates of creation in the case of some image
% files using exif info).
% All optional input parameters must be specified as parameter/value pairs,
% e.g. as in 
%          writeFileList2TextFile(dDir,outFn,'fileExt',{'jpg'})
%
%                         >>> INPUT VARIABLES >>>
% NAME         TYPE/DEFAULT        DESCRIPTION
% dDir         char                directory containing files
% outFn        char                path & name of output text file
% fileExt      cell, {'*'}         extensions of files to be listed, e.g.
%                                  {'jpg'}
% doSubdir     logical, true       if true, subdirs will be searched
% sortCrit     char, 'none'        'none' - files will not be sorted
%                                  'date' - date of creation/alteration
%                                  'dateEXIF' - date of creation
%                                   according to EXIF data (image files)
%
%                         <<< OUTPUT VARIABLES <<<
% NAME           TYPE/DEFAULT           DESCRIPTION
%
%

%                 - work in progress - 

% to do
% - make this into a much more flexible and powerful function which collects
% also interesting meta-info like the device used to take the picture

sortCrit='none';
fileExt={'*'};
doSubdir=true;
pvpmod(varargin,{'fileExt','doSubdir','sortCrit'})

% subdirectories
if doSubdir
  masterFileList=dir(dDir);
  subdirList=masterFileList([masterFileList.isdir]);
  subdirList=setdiff({subdirList.name},{'..'});
else
  subdirList={'.'};
end

fileList=dir([dDir subdirList{1} '/*.' fileExt{1}]);
for g=2:numel(subdirList)
  for k=1:numel(fileExt)
    fileList=cat(1,fileList,dir([dDir subdirList{g} '/*.' fileExt{k}]));
  end
end
numImg=numel(fileList);

switch sortCrit
  case 'dateEXIF'
    imgCreationDate=NaT(numImg,1);
    % retrieve exif info
    for k=1:numel(fileList)
      imgInfo=imfinfo([fileList(k).folder filesep fileList(k).name]);
      if isfield(imgInfo,'DateTime')
        % include fallbacks if specific format is not met, and cushion in
        % try-catch
        imgCreationDate(k)=datetime(imgInfo.DateTime,'InputFormat','yyyy:MM:dd HH:mm:ss');
      else
        warning(['file ' fileList(k).name ' has no creation date (likely from Whatsapp) - using modification date']);
        % imgCreationDate(k)=datetime(imgInfo.FileModDate,'InputFormat','dd-mmm-yyyy HH:MM:SS');
        imgCreationDate(k)=datetime(imgInfo.FileModDate);
      end
    end
    % sort
    [~,sortIx]=sort(imgCreationDate);

  case 'none'
    sortIx=1:numImg;
  
  otherwise
    error('bad sortCrit')
end
fileList=fileList(sortIx);

% write list
fid=fopen(outFn,'w');
for k=1:numel(fileList)
  fprintf(fid,'%s\n',[fileList(k).folder filesep fileList(k).name]);
end
fclose(fid);