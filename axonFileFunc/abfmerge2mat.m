function [mergeFn,fi,varargout]=abfmerge2mat(fn,varargin)
% ** function [mergeFn,fi,varargout]=abfmerge2mat(fn,varargin)
% merges abf data files, heeding timing, optionally treats them in multiple
% ways (filtering, spike removal, base line subtraction) and saves data as
% singles to matfile(s). Each channel in the abf file will be represented
% by a separate variable of the same name excluding blanks (e.g. 'IN 0'
% will result in variable IN0). 
% The function was primarily intended to produce downsampled LFP traces of
% multiple abf files, but can also produce "MUA" data, that is, a signal
% reflecting (multi-unit) spiking activity, provided proper filter
% frequencies are specified in input variable cFreq. For this purpose, the
% data will first be bandpass filtered, then the absolute value is taken,
% and finally a lowpass filter is applied. This mode of signal treatment is
% engaged if optional input argument 'muaCFreq' is nonempty.
% !! It is assumed that neither channels nor sampling interval nor any 
% other viable variable changes from one recording to the next !!
% 
%                         >>> INPUT VARIABLES >>>
% NAME        TYPE,DEFAULT             DESCRIPTION
% fn          cell arr                 file names WITHOUT extension
% dDir        char arr, ''             directory in which files reside    
% ch          cell arr or 'all','all'  channels to be read & merged
% cFreq       array, [nan nan]         edge frequencies (Hz) of high- and
%                                       lowpass filters
% muaCFreq    scalar,[]                if nonempty, "mua" operating mode is
%                                       engaged (see above) and the value
%                                       is the edge frequency of the final
%                                       lowpass filter
% lineFreq    array, nan               center frequency (Hz) of line pickup
%                                        to be eliminated (multiple 
%                                        frequencies allowed)
% fiFunc      char arr, 'filtfilt'     filter function to use: either
%                                        'filtfilt' or 'filter'
% sampFreq    scalar, nan              the sampling frequency (Hz) to which 
%                                       data will be up- or downsampled
% doBaseSub   array, []                base line (mean of all points) will 
%                                        be subtracted from all channels 
%                                        represented by a 1
% gapVal      scalar, 0                value to which data points in gaps 
%                                        between recordings will be set
% noMerge     scalar, 0                if nonzero files will not be merged
%                                        and no data will be returned but 
%                                        the name of the resulting file as
%                                        if the files were merged will be 
%                                        returned 
% tslC        cell arr of tls          list of events (e.g. spx) to be
%                                        eliminated from raw data prior to
%                                        filtering. Channels in columns,
%                                        files in rows.
% sIntvC      cell arr of arr, []      corresponding intervals in which to 
%                                        'eliminate' spikes 
% fnAdd       char, ''                 addition to name of resulting mat
%                                        file
% 
%                         <<< OUTPUT VARIABLES <<<
% NAME        TYPE/DEFAULT             DESCRIPTION
% mergeFn                              name of the merged file
% fi                                   file information (** NOTE:
%                                      channel-specific information is
%                                      given only for requested channels
%                                      but in the order inherent to the abf
%                                      file, see below)
% varargout                            the data written into file
%                                      (memory-consuming!) ** NOTE:
%                                      channel order ignores the order in
%                                      'ch', and instead always matches the 
%                                      order inherent to the abf file, to 
%                                      be retrieved in output variable fi

% §§ to do: 
% - think about fields of fi worth saving
% - sampFreq is nan -> use original freq
% - it should be possible to set filtering & base line subtraction channel
% by channel

% --- defaults
dDir='';
ch='all';
cFreq=[nan nan];
muaCFreq=[];
lineFreq=nan;
fiFunc='filtfilt';
doBaseSub=[];
sampFreq=nan;
gapVal=0;
noMerge=0;
tslC=[];
sIntvC=[];
fnAdd='';
pvpmod(varargin)

% --- check input & create handy parameters
if ~isempty(dDir)
  dDir=[dDir '\'];
end
switch fiFunc
  case 'filtfilt'
    force_sigproc=1;
  case 'filter'
    force_sigproc=0;
  otherwise
    error('bad choice for input parameter fiFunc');
end

stop='e';
% if only a single file was specified as a char convert it to a cell
if ischar(fn)
  fn={fn};
end
if numel(unique(fn))<numel(fn)
  error('duplicate files');
end
nF=length(fn);
if nF>1
  mergeFn=[fn{1} '_'  fn{end}(end-3:end) fnAdd '.mat'];
else
  mergeFn=[fn{1} fnAdd '.mat'];
end

if ~isempty(muaCFreq) && isfinite(muaCFreq)
  if ~all(isfinite(cFreq))
    error('cFreq must be fully specified for MUA processing');
  end
end

if ~isempty(tslC)
  doEliminateSpx=true;
  if ~iscell(tslC)
    % most likely user forgot to wrap single tsl into a cell
    warning('tslC is supposed to be a cell array');
    tslC={tslC};
  end
  [n1,n2]=size(tslC);
  if n1~=nF
    error('number of files and number of rows in input var tslC must match');
  end
  % correct number of channels = columns in tslC cannot yet be checked
else
  doEliminateSpx=false;
end

if ~isempty(sIntvC)
  if ~doEliminateSpx
    error('sIntvC was specified as input, but tslC is missing')
  end
  if ~iscell(sIntvC)
    % most likely user forgot to wrap single sIntvC into a cell
    warning('sIntvC is supposed to be a cell array');
    sIntvC={sIntvC};
  end
else
  if doEliminateSpx;
    error('tslC was specified as input, but sIntvC is missing')
  end
end

  
if nargout==3
  varargout{1}=[];
  doReturnD=true;
else
  doReturnD=false;  
end

% determine total length of merged streams
[nix,nix2,fi1]=abfload([dDir fn{1} '.abf'],'info');
[nix,nix2,fi2]=abfload([dDir fn{end} '.abf'],'info');

% channels
if ischar(ch) 
  if ~strcmpi(ch,'all')
    error('ch must be a cell array containing channel names or string ''all''');
  else
    ch=fi1.recChNames;
    chIx=1:numel(ch);
    tmpIx=chIx;
  end
else
  [nix,chIx,tmpIx]=intersect(fi1.recChNames,ch);
  % ** the following two lines ensure that channel order is according to
  % the internal abf order as listed in the header
  chIx=sort(chIx);
  ch=fi1.recChNames(chIx);
  % furthermore, rearrange tslC, sIntvC and doBaseSub accordingly
  if ~isempty(tslC)
    tslC=tslC(:,tmpIx);
  end
  if ~isempty(sIntvC)
    sIntvC=sIntvC(:,tmpIx);
  end
end

% now check & expand doBaseSub
if isempty(doBaseSub)
  doBaseSub=nan(size(ch));
elseif numel(doBaseSub)==1
  doBaseSub=repmat(doBaseSub,size(ch));
elseif numel(doBaseSub) ==numel(ch)
  doBaseSub=doBaseSub(tmpIx);
else
  error('input argument ''doBaseSub'' must contain one value or as many values as there are channels in the abf file');
end


% fi is a 'made-up' structure containing only vital information about the
% merged file based on header information in the first file. 
fi.origSi=fi1.si;
if isfinite(sampFreq)
  fi.si=1e6/sampFreq;
else
  fi.si=fi1.si;
end
% fi.recTime is specified in seconds
fi.recTime=[fi1.recTime(1) fi2.recTime(2)];
fi.OrigDataPtsPerChan=cont2discrete(fi2.recTime(2)-fi1.recTime(1),fi.origSi*1e-6,'intv',1);
fi.dataPtsPerChan=cont2discrete(fi2.recTime(2)-fi1.recTime(1),fi.si*1e-6,'intv',1);
fi.lFileStartTime=fi1.lFileStartTime;
% put out only those channels that are actually chosen
fi.recChUnits=fi1.recChUnits(chIx);
fi.recChNames=fi1.recChNames(chIx);

% return if user requests only information on merged files 
if noMerge
  return
end

% are the files sweep-based?
isSweepType=false;
% preallocate
d=repmat(gapVal,fi.dataPtsPerChan,1);

% do things channel by channel, then file by file (less memory demand than file by file)
for chInd=1:length(ch)
  dbch=ch{chInd};
  disp(['channel: ' dbch '...']);
  % deblank names
  dbch=dbch(~isspace(dbch));
  for fIx=nF:-1:1
    [nix,nix2,fiCurr]=abfload([dDir fn{fIx} '.abf'],'info');
    [tmpD,si]=abfload([dDir fn{fIx} '.abf'],'start',0,'stop',stop,'channels',ch(chInd));
    % if tmpD is a 3D array we're dealing with a recording in fixed-length
    % mode which is supposed to be concatenated. Using information from
    % abfload glue as many gapVals to the bottom of the sweeps as are
    % necessary to represent the gaps between sweeps
    if ismember(fiCurr.nOperationMode,[2 5])
      if fIx~=nF && ~isSweepType
        error('files with differing recording modes may not be concatenated')
        % § in fact, that may be doable, but requires close scrutiny
      end
      isSweepType=true;
      if doEliminateSpx
        error('the present version does not support spike elimination in sweep-based recording modes');
        % § again, this is implementable, see e.g. threshdetgui
      end
      % due to an imperfection in abfload these values may be noninteger
      fiCurr.sweepStartInPts=round(fiCurr.sweepStartInPts);
      tmpStartDiff=unique(diff(fiCurr.sweepStartInPts));
      if numel(fiCurr.sweepStartInPts)>1
        if numel(tmpStartDiff)==1
          tmpD=cat(1,tmpD,repmat(gapVal,[tmpStartDiff-fiCurr.sweepLengthInPts  1  fiCurr.lActualEpisodes]));
          % now make into 1D column array
          tmpD=tmpD(:);
        else
          error('cannot concatenate sweeps because start times are not regularly spaced')
        end
      end
    end
    
    % kick out events if requested
    if doEliminateSpx
      sIntv=sIntvC{fIx,chInd};
      if all(isfinite(sIntv))
        % default for interpolation: substitution interval flanked by 2 ms
        % on each side
        eIntv=sIntv+[-2 2];
        tsl=tslC{fIx,chInd};
        if ~isempty(tsl) && any(isfinite(tsl))
          tmpD=etslexcsubst(tmpD,si,tsl,sIntv,eIntv);
        else
          disp('empty tsl, no spike substitution');
        end
      else
        disp('sIntv set to nan, no spike substitution');
      end
    end
    
    if ~isempty(muaCFreq) && isfinite(muaCFreq)
      % 'mua' mode: 
      % - bandpass
      tmpD=bafi(tmpD,si,cFreq,'force_sigproc',force_sigproc);
      % - abs and lowpass
      tmpD=lofi(abs(tmpD),fi.si,muaCFreq);
      % - resample to target freq
      if isfinite(sampFreq)
        nPt=size(tmpD,1);
        tmpD=interp1((1:nPt)',tmpD,linspace(1,nPt,nPt*si*(sampFreq/1e6))');
      end
    else
      % lowpass filter
      if isfinite(cFreq(2))
        tmpD=lofi(tmpD,si,cFreq(2),'force_sigproc',force_sigproc);
      end
      % resample to target freq
      if isfinite(sampFreq)
        nPt=size(tmpD,1);
        tmpD=interp1((1:nPt)',tmpD,linspace(1,nPt,nPt*si*(sampFreq/1e6))');
      end
      % highpass filter
      if isfinite(cFreq(1))
        tmpD=hifi(tmpD,fi.si,cFreq(1),'force_sigproc',force_sigproc);
      end
      % line pickup elimination
      if all(isfinite(lineFreq))
        tmpD=elim_hum(tmpD,fi.si,lineFreq);
      end
    end
    
    % base line subtraction
    if doBaseSub(chInd)==1
      tmpD=tmpD-mean(tmpD);
    end
    % allocate right slots..
    ix=cont2discrete(fiCurr.recTime-fi1.recTime(1),fi.si*1e-6,'intv',1);
    switch diff(ix)+1-length(tmpD)
      case 0
        ;
      case -1
        disp('real number of data points one more than expected');
        tmpD(end)=[];
      case 1
        disp('real number of data points one less than expected');
        tmpD(end+1)=tmpD(end);
      otherwise
        error('sampling time index not OK');
    end
    % embed
    d(ix(1):ix(2))=tmpD;
  end
  % save data as singles in variable with deblanked channel name
  eval([dbch '=single(d);']);
  % save in matfile
  if chInd==1,
    save([dDir mergeFn],dbch,'fi','-mat');
  else
    save([dDir mergeFn],dbch,'-mat','-append');
  end
  % delete
  eval(['clear ' dbch ';']);
  % collect data if requested
  if doReturnD
    varargout{1}=[varargout{1} d];
  end
end



