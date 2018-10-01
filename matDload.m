function [dat,si,fi]=matDload(fn,varargin)
% ** function [dat,si,fi]=matDload(fn,varargin)
% loads time series data from *.mat files, mimicking the behavior of
% abfload with gap-free data.
% The time series data has to be in a specific format and some additional
% information on the data must reside in the file.
% RAW DATA:
% Each channel must reside in a separate double array. No restriction on 
% variable name(s) except that 'dat' is forbidden.
% OTHER INFORMATION 
% The *.mat file must contain a struct array named either 'abfi' (historic
% reasons) or 'fi' (for 'file information').
% This struct array must contain the following fields:
% - si (the sampling interval in microseconds)
% - dataPtsPerChan (number of data points per channel)
% Also often needed is 
% - recTime (recording start and stop times in seconds from midnight)
%
% All optional input parameters listed below (= all except the file name) 
% must be specified as parameter/value pairs, e.g. as in 
%          dat=matDload('d:\data01.mat','channels',{'ch1','ch2'},'stop',70);
%
%                    >>> INPUT VARIABLES >>>
%
% NAME        TYPE, DEFAULT      DESCRIPTION
% fn          char array         data file name
% channels    cell array         names of channels to be read, like {'ch0','ch8'};
% start       scalar, 0          start of excerpt to be read (unit: sec)
% stop        scalar or char,    end of excerpt to be read (unit: sec). 
%             'e'                 May be set to 'e' (end of file).
%
%                         <<< OUTPUT VARIABLES <<<
% NAME        TYPE/DEFAULT       DESCRIPTION
% dat                            the data written into file
%                                ** NOTE: channel order ignores the order
%                                in 'channels', and instead always matches
%                                the order inherent to the abf file, to be
%                                retrieved in output variable fi
% si                             sampling interval in us
% fi                             file information


% defaults   
start=0.0;
stop='e';
channels={''};
pvpmod(varargin);

if strcmpi(channels,'a'),
  error('loading data from matfile: channel names must be given explicitly');
elseif isempty(channels{1})
  error('input variable ''channels'' must be specified');
end

% load file information (abfi in older versions, fi in newer versions)
s=whos('-file',fn);
varName={s.name};
if any(strcmp('fi',varName))
  load(fn,'fi');
  if any(strcmp('abfi',varName))
    error('data file must contain either ''abfi'' or ''fi'', but not both');
  end
elseif any(strcmp('abfi',varName))
  load(fn,'abfi');
  fi=abfi;
else
  error('data file must contain either ''abfi'' or ''fi''');
end

% --- deblank all channel names
deblankRecChNames=fi.recChNames;
for cIx=1:numel(deblankRecChNames)
  deblankRecChNames{cIx}=deblankRecChNames{cIx}(~isspace(deblankRecChNames{cIx}));
end
for cIx=1:numel(channels)
  channels{cIx}=channels{cIx}(~isspace(channels{cIx}));
end

% --- sort requested channels according to order defined in header
tmpIx=nan(numel(channels),1);
for ci=1:numel(channels)
  tmpIx(ci)=find(strcmp(channels(ci),deblankRecChNames));
end
[nix,tmpIx]=sort(tmpIx);
channels=channels(tmpIx);

% --- deal with start & stop times
ix(1)=cont2discrete(start*1e6,fi.si,'intv',0);
if strcmpi(stop,'e')
  ix(2)=fi.dataPtsPerChan;
else
  ix(2)=cont2discrete(stop*1e6,fi.si,'intv',1);
end
% ---- load
% preallocate
dat=repmat(nan,diff(ix)+1,length(channels));
for chInd=1:length(channels)
  dbch=channels{chInd};
  load(fn,dbch);
  eval(['dat(:,chInd)=double(' dbch '(ix(1):ix(2),:));']);
  eval(['clear ' dbch ';']);
end
si=fi.si;
