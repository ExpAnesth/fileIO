function [d,si,h]=read_mcdhdf5(fn,varargin)
% ** function [d,si,h]=read_mcdhdf5(fn,varargin)
% returns AnalogStream data from hdf5-converted mcd (MCRack data) files as
% regular double or single arrays, making use of the McsHDF5 Toolbox. It is
% modeled on abfload.m. If the second input variable is the char array
% 'info' as in
%         [d,si,h]=read_mcdhdf5('d:\data01.h5','info') 
% the function will not read any AnalogStream data but only return detailed
% information on the file in output variable h; d and si will be empty. In
% all other cases read_mcdhdf5 will read and return AnalogStream data.
% Optional input parameters listed below (that is, all except fn) must be
% specified as parameter/value pairs, e.g. as in
%         d=read_mcdhdf5('d:\data01.h5','start',100,'stop','e');
% **NOTE: extracting the data from McsHDF5.xxx objects results in spikes of
% memory demand which may be substantial, up to twice the size of the total
% analog data in the file, depending on the reading method. If you
% encounter out of memory errors with large data files, set optional input
% parameter doReadSeq to true (reading will be excruciatingly inefficient
% and slow, but this may be the only way to retrieve and hold the requested
% data in memory as ). Alternatively, read data in smaller blocks of contiguous
% channel groups and concatenate.
%
% >>> INPUT VARIABLES >>>
% NAME        TYPE, DEFAULT      DESCRIPTION
% fn          char array         hdf5 data file name
% start       scalar, 0          start of cutout to be read (unit: s)
% stop        scalar or char,    end of cutout to be read (unit: s). May be
%             'e'                 set to 'e' (end of file).
% channels    cell array         names of channels to be read, like 
%              or char, 'a'       {'32','64'}. Make sure spelling is
%                                 correct, including blanks. If set to 'a',
%                                 all channels will be read.
%                                 *****************************************
%                                 NOTE: channel order in output variable d
%                                 ignores the order in 'channels', and
%                                 instead always matches the order inherent
%                                 to the hdf file, to be retrieved in
%                                 output variables h and mcsd!
%                                 *****************************************
% recIx       scalar, 1          index into .Recording
% anStreamIx  scalar, 1          index into .Recording{recIx}.AnalogStream
% dataType    char, 'double'     'double' or 'single'; the latter may be
%                                 advisable if requested data chunk
%                                 approaches memory available to Matlab
% signalExp   char, 'm'          desired exponent of signal unit, options
%                                 are '','m','µ'
% doReadSeq   logical, false     if true, slow but memory-saving sequential 
%                                 method of reading data is chosen even if
%                                 (almost) all data is requested
% doDispInfo  logical, true      if true, information on the loaded file
%                                 will be put out to console (if false,
%                                 only ínformation on erroneous input will
%                                 be displayed)
% << OUTPUT VARIABLES <<<
% NAME  TYPE            DESCRIPTION
% d                     the data read, <data pts> by <number of chans>
% si    scalar          the sampling interval in µs
% h     struct          information on file (selected header parameters 
%                        modeled on those produced by abfload.m)

% -------------------------------------------------------------------------
% PART 1: check of input vars
% -------------------------------------------------------------------------
% --- defaults   
start=0.0;
stop='e';
channels='a';
recIx=1;
anStreamIx=1;
dataType='double';
signalExp='m';
doReadSeq=false;
doDispInfo=true;
% if first and only optional input argument is string 'info' the user's
% request is to obtain information on the file (header parameters), so set
% flag accordingly
if nargin==2 && ischar(varargin{1}) && strcmp('info',varargin{1})
  doLoadData=false;
  % if no output argument is requested assume that the user would like to
  % obtain some basic information on the file on the command window, so
  % leave doDispInfo at its default (true value). Do the same if only one
  % output arg is specified (which is most certainly done inadvertently
  % because in 'info' mode this will be an empty array). In all other cases
  % assume that text output is not required so suppress it
  if nargout>1
    doDispInfo=false;
  end
else
  doLoadData=true;
  % assign values of optional input parameters if any were given
  pvpmod(varargin);
end

% check dataType and, doing so, set up cfg
cfg=[];
if isempty(intersect(dataType,{'single','double'}))
  error('illegal dataType')
else
  cfg.dataType=dataType;
end

% output variables
d=[]; 
si=[];
h=[];

% -------------------------------------------------------------------------
% PART 2a: create header parameters
% -------------------------------------------------------------------------
% ** Note: the parameters created here are modeled after the ones produced
% by function abfload.m, so some of them make sense only viewed in this
% context
dispif(doDispInfo,['opening ' fn '...']);
mcsd=McsHDF5.McsData(fn,cfg);
% sampling interval
si=double(unique(mcsd.Recording{recIx}.AnalogStream{anStreamIx}.Info.Tick));
if numel(si)>1
  error('nonuniform sampling intervals across channels');
end
h.si=si;
% recording start time in seconds from midnight:
% fields .TimeStamp and .DateInTicks are dates given as the number of
% 100-nanosecond intervals since January 1, 0001, so they need to be
% converted with ticks2days and datetime
dateInDays=datetime(ticks2days(mcsd.Recording{recIx}.TimeStamp),'ConvertFrom','datenum');
% compute the number of seconds passed since midnight
h.lFileStartTime=seconds(dateInDays-dateshift(dateInDays,'start','day'));
% start and stop time of recording in seconds after midnight (** note that
% .Duration is given in microseconds, yeah)
h.recTime=h.lFileStartTime+[0 double(mcsd.Recording{recIx}.Duration)/1e6];
% names of channels
h.recChNames=mcsd.Recording{recIx}.AnalogStream{anStreamIx}.Info.Label;
% number of recorded channels
h.nADCNumChannels=numel(h.recChNames);
% physical units of recorded signal (e.g. mV), composed of user-defined
% exponential and unit defined by recording
h.recChUnits=cellfun(@horzcat,repmat({signalExp},h.nADCNumChannels,1),...
  mcsd.Recording{recIx}.AnalogStream{anStreamIx}.Info.Unit,...
  'Uniformoutput',false);
% set mFactor, by which signals have to be multiplied to result in the
% user-defined unit
switch signalExp
  case ''
    mFactor=1;
  case 'm'
    mFactor=1e-3;
  case 'µ'
    mFactor=1e-6;
  otherwise
    error('bad signalExp')
end
mFactor=10.^double(mcsd.Recording{recIx}.AnalogStream{anStreamIx}.Info.Exponent)./mFactor;
% number of data points per channel
h.dataPtsPerChan=numel(mcsd.Recording{recIx}.AnalogStream{anStreamIx}.ChannelDataTimeStamps);
% total number of data points in file
h.dataPts=h.dataPtsPerChan*h.nADCNumChannels;

% -------------------------------------------------------------------------
%  PART 2b: compute derived parameters & perform some plausibility checks
% -------------------------------------------------------------------------
% deal with time excerpt
if ischar(stop)
  if ~strcmpi(stop,'e')
    error('input parameter ''stop'' must be specified as ''e'' (=end of recording) or as a scalar');
  else
    % .Duration is given in microseconds
    stop=double(mcsd.Recording{recIx}.Duration)/1e6;
  end
end

% set up flag deciding which data reading method to use: if user explicitly
% requested sequential reading of channels, or if requested time interval
% is <= 25% of full recording interval use partial reading routine as it's
% likely faster and definitely memory-saving compared to a full read and
% subsequent pruning
doReadPartial= doReadSeq || (stop-start)/(double(mcsd.Recording{recIx}.Duration)/1e6)<=.25;

% the numerical value of all recorded channels (numbers 0..60 or so,
% depending on the recording hardware)
recChIdx=mcsd.Recording{recIx}.AnalogStream{anStreamIx}.Info.ChannelID;
% the corresponding indices into loaded data d
recChInd=1:length(recChIdx);

% check whether requested channels exist
eflag=0;
if ischar(channels)
  if strcmp(channels,'a')
    chInd=recChInd;
  else
    fclose(fid);
    error('input parameter ''channels'' must either be a cell array holding channel names or the single character ''a'' (=all channels)');
  end
else
  % check for requested channels which do not exist
  missingChan=setdiff(channels,h.recChNames);
  % identify requested channels among available ones
  [~,chInd]=intersect(h.recChNames,channels);
  % ** index chInd must be sorted because intersect sorts h.recChNames
  % alphanumerically, which needs not necessarily correspond to the order
  % inherent in the hdf file
  chInd=sort(chInd);
  if isempty(chInd) || ~isempty(missingChan)
    % set error flag to 1
    eflag=1;
  end
end
if eflag
  disp('**** available channels:');
  disp(h.recChNames);
  disp(' ');
  disp('**** requested channels:');
  disp(channels);
  error('at least one of the requested channels does not exist in data file (see above)');
else
  % the names of the returned channels
  h.readRecChNames=h.recChNames(chInd);
  % the number of read channels
  nReadChan=numel(chInd);
  dispif(doDispInfo,'**** reading channels:');
  dispif(doDispInfo,h.readRecChNames);
end

% modify flag deciding which data reading method to use: if requested
% number of channels is less than 90% of available channels use partial
% reading routine as it's likely faster and definitely memory-saving
% compared to a full read and subsequent pruning
doReadPartial=doReadPartial || nReadChan/h.nADCNumChannels<=.9;

% -------------------------------------------------------------------------
%    PART 3: read data 
% -------------------------------------------------------------------------
dispif(doDispInfo,['total length of recording: ' num2str(diff(h.recTime),'%5.1f') ' s ~ ' num2str(diff(h.recTime)/60,'%3.0f') ' min']);
dispif(doDispInfo,['sampling interval: ' num2str(h.si,'%5.0f') ' µs']);
if doLoadData
  if doReadPartial
    % if readPartial functionality is used, define cfg.window so the data
    % won't have to be curbed afterwards
    cfg.window=[start stop];
    % if a contiguous block of channels is requested and sequential reading
    % is not explicitly requested...
    if (nReadChan==1 || all(diff(chInd)==1)) && ~doReadSeq
      cfg.channel=chInd([1 end])';
      % note that this overwrites the metadata in mcsd
      mcsd=mcsd.Recording{recIx}.AnalogStream{anStreamIx}.readPartialChannelData(cfg);
      % ** note that
      % 1. assigning the values of field .ChannelData from a
      % McsHDF5.McsAnalogStream object to an array is very awkward because
      % we will have two copies of the channel data in memory until the
      % object is deleted, but there is no alternative if output variable d
      % is to be an array
      % 2. the code makes use of the 'new' automatic expansion of arrays
      % with elementwise multiplication
      % 3. data need to be transposed
      d=(mcsd.ChannelData.*mFactor(chInd))';
      clear mcsd
    else
      % no alternative to a loop here
      for k=1:nReadChan
        cfg.channel=chInd([k k])';
        tmpMcsd=mcsd.Recording{recIx}.AnalogStream{anStreamIx}.readPartialChannelData(cfg);
        if k==1
          % make sure that d is of correct class via non-indexed assignment
          d=tmpMcsd.ChannelData'*mFactor(chInd(k));
          numDataPoints=numel(d);
          % preallocate via automatic expansion
          d(numDataPoints,nReadChan)=0;
        else
          d(1:numDataPoints,k)=tmpMcsd.ChannelData'*mFactor(chInd(k));
        end
      end
      clear tmpMcsd mcsd
    end
  else
    % compute time indexes (**note: readPartialChannelData puts out more
    % data points than seems reasonable, but do not heed this and instead
    % compute the indexes properly here)
    winIx=cont2discrete([start stop],si/1e6,'intv',true);
    % ** note that this is the point of peak memory demand
    d=mcsd.Recording{recIx}.AnalogStream{anStreamIx}.ChannelData;
    clear mcsd
    % d is channels x time!
    d=d(chInd,winIx(1):winIx(2));
    d=(d.*mFactor(chInd))';
  end  
end

% ============================ LOCAL FUNCTIONS ============================

function dispif(doDispInfo, msg)
if doDispInfo
  disp(msg)
end