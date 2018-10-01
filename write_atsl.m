function write_atsl(atsl,fn,varargin)
% ** function write_atsl(atsl,fn,varargin)
% writes an 'advanced time stamp list' into a text file.
% All HEADER input parameters listed below are optional and must be 
% specified as parameter/value pairs, e.g.  
%          write_atsl(tsl,'2007_10_02_0001_IN1.txt','rectime',180);
%
%                         >>> INPUT VARIABLES >>>
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% atsl             2-column array        see define_atsl.m
% fn               char array            name of file in which atsl shall
%                                          be written (including extension)
% varargin         parameter/value pairs as listed below:
% *************************************************************************
% **   THE FOLLOWING ARE HEADER PARAMETERS AS DEFINED IN DEFINE_ATSL     **
% *************************************************************************
% atsl_version
% rectime
% ephys_filename
% ephys_chan
% ev_threshold
% burst_interval1
% burst_interval2
% burst_minev
% silentper_dur
% silentper_maxev
% data_quality
% comment

% ----- default values of header parameters
% the following lines of code generate new variables, namely all header
% parameters by which an atsl is defined, and assign default values to them
atslDef=define_atsl;
nField=size(atslDef,1);
for g=1:nField
  eval([atslDef{g,1} '=atslDef{g,2};']);
end

% ----- checks of input of user
if nargin<2
  error([mfilename ' requires at least two input arguments']);
end
% - size of atsl
[n1,n2]=size(atsl);
if n1==0
  warning('atsl is empty');
  isEmptyAtsl=true;
else
  isEmptyAtsl=false;
end
if ~isEmptyAtsl && n2~=2
  error('atsl must contain two columns');
end
% - time stamps must be sorted
if ~isEmptyAtsl && any(diff(atsl(:,1)))<=0
  error('time stamps in atsl are not sorted');
end

% format string for writing time stamps and burst information: for values
% of -1 and 0 in the second column, use no decimal digits, otherwise use
% four for both time stamps and burst info
tsFormatString_nodecimal='%9.4f %8.0f\n';
tsFormatString='%9.4f %8.4f\n';

% - fn must be a nonempty char arr
if isempty(fn) || ~ischar(fn)
  error('check fn');
end

% - if any of the header parameters are misspelled an error will follow
badSheep=setdiff(varargin(1:2:end),atslDef(:,1));
if ~isempty(badSheep)
  disp(' ')
  disp(badSheep);
  disp(' ')
  error('the input parameter(s) listed above do not do not exist in the atsl V. 1.0 standard');
end

% now we can assign values specified by the user to the variables
pvpmod(varargin);

% do the inverse of the above - collect values of parameters in atslDef so
% printing them in the output file can be done without explicitly referring
% to their name
for g=1:nField
  eval(['atslDef{g,2}=' atslDef{g,1}  ';']);
end

fid=fopen(fn,'wt');
% write the header
for g=1:nField
  fprintf(fid,'%% %15s\t', atslDef{g,1});
  fprintf(fid,[atslDef{g,3} '\n'],atslDef{g,2});
end

% if atsl is empty the job is done here
if ~isEmptyAtsl
  try
    % §§§ blockwise writing of atsl has yet to undergo thorough check
    % index to lines whose second column shall be written with decimal
    % digits
    ix=find(atsl(:,2)~=0 & atsl(:,2)~=-1);
    % write atsl in blocks
    % ** note that we have to transpose atsl before writing **
    % - ts before first burst
    if ix(1)>1
      fprintf(fid,tsFormatString_nodecimal,atsl(1:ix(1)-1,:)');
    end
    % - all following
    for g=1:numel(ix)-1
      % - single burst line
      fprintf(fid,tsFormatString,atsl(ix(g),:)');
      % - following ts, if any
      if ix(g)<n1
        fprintf(fid,tsFormatString_nodecimal,atsl(ix(g)+1:ix(g+1)-1,:)');
      end
    end
    % - last burst line
    fprintf(fid,tsFormatString,atsl(ix(end),:)');
    % - following ts, if any
    if ix(end)<n1
      fprintf(fid,tsFormatString_nodecimal,atsl(ix(end)+1:n1,:)');
    end
  catch
    warndlg('writing burst information of atsl in proper format failed - now writing in simple format (there will be no loss of information)');
    fprintf(fid,tsFormatString,atsl');
  end
end
fclose(fid);