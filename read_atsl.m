function [head,atsl,varargout]=read_atsl(fn)
% ** function [head,atsl,varargout]=read_atsl(fn)
% reads an 'advanced time stamp list' (atsl) from a text file and puts out
% the header parameters and atsl as well as time stamp lists derived
% from it.
% Optional output parameters must be explicitly specified as in e.g.  
%          [head,atsl,etsl]=read_atsl('2007_10_02_0001_IN1.txt');

%                         >>> INPUT VARIABLES >>>
% NAME             TYPE/DEFAULT          DESCRIPTION
% fn               char                  atsl file name 
%
%                         >>> OUTPUT VARIABLES >>>
% NAME             TYPE/DEFAULT          DESCRIPTION
% head             struct                header parameters - see define_atsl.m 
% atsl             2 col array           atsl
% varargout{1}     2D array              etsl
% varargout{2}     1 col array           tsl

% start by informing us which file we're dealing with
disp(['** reading ' fn '...']);

isEmptyAtsl=false;
etslconst;

% retrieve info about definition of atsl
atslDef=define_atsl;
nField=size(atslDef,1);
% set up header var with default values as defined in define_atsl
head=cell2struct(atslDef(:,2),atslDef(:,1),1);

% general strategy for reading data:
% - read header (MUST be nField lines)
% - check existence of each parameter and assign value
% - read num data, if any exist
% - if so, check num data

fid=fopen(fn,'r');
% simplistic header reading & integrity check:
% - read nField lines
% - check whether all begin with a '%' and whether all parameters listed in
% define_atsl exist 
% - if that is not so consider header corrupt and issue an error
% - if everything is fine read numbers
% - if fscanf fails issue error
parIx=[];
for g=1:nField
  % get line...
  s=fgetl(fid);
  % check whether the first element is the percentage sign
  if ~strcmp(s(1),'%')
    fclose(fid);
    error(['header is corrupt (line ' int2str(g) ' should begin with a ''%''']);
  end
  % read name of parameter (which by default contains no whitespace) and
  % its value
  parNm=sscanf(s(2:end),'%s',1);
  % check whether parameter name is legal and find its index into atslDef
  tmpIx=strmatch(parNm,atslDef(:,1),'exact');
  if isempty(tmpIx)
    fclose(fid);
    error(['header contains illegal parameter ''' parNm ''])
  else
    parIx=[parIx tmpIx];
  end
  % the format string in atslDef (3rd column) can be used for determining
  % whether header parameter is numeric or char...
  if any(ismember(atslDef{tmpIx,3}(end),{'s','c'}))
    parVal=fliplr(deblank(fliplr(s(strfind(s,parNm)+length(parNm):end))));
  elseif any(ismember(atslDef{tmpIx,3}(end),{'e','f','g'}))
    parVal=sscanf(s(2:end),'%*s %f');    
  else
    fclose(fid);
    error('define_atsl contains an illegal format string - tell the programmer...');
  end
  head=setfield(head,parNm,parVal);
end
% final check: all header parameters present?
sDiff=setdiff(1:nField,parIx);
if ~isempty(sDiff)
  fclose(fid);
  errordlg({'header parameters are missing: ', strvcat(atslDef{sDiff,1})});
  error('see error window')
end

% finally, read time stamps
[atsl,nTs]=fscanf(fid,'%f %f',[2,inf]);
if nTs==0
  [msg,errnum]=ferror(fid);
  fclose(fid);
  if errnum==-4
    warning('atsl contains no events');
    % make sure atsl contains proper number of columns
    atsl=zeros(0,2);
    isEmptyAtsl=true;
  else
    error(['reading of time stamps failed - here is the error message by ferror:' msg]);
  end
else
  fclose(fid);
  % don't forget to flip
  atsl=atsl';
end

if isEmptyAtsl
  if nargout>2
    varargout{1}=zeros([0 max([etslc.tsCol etslc.durCol])]);
  end
  if nargout>3
    varargout{2}=zeros([0 1]);
  end
else
  % check time stamp list
  [n1,n2]=size(atsl);
  if n2~=2
    error('atsl must contain two columns');
  end
  % time stamps must be sorted
  if any(diff(atsl(:,1)))<=0
    error('time stamps in atsl are not sorted - possibly header is corrupt');
  end
  % finally, split up, if so desired
  if nargout>2
    tmpIx=logical(atsl(:,2));
    % ** etsl and tsl have a time basis of ms, so convert here
    varargout{1}(1:length(find(tmpIx)),[etslc.tsCol etslc.durCol])=atsl(tmpIx,:)*1000;
  end
  if nargout>3
    varargout{2}=atsl(:,1)*1000;
  end
end