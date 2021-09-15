function [ ids, times ] = load_netsim_spikes( filepath )
% *NETSTIM*
%
%	LOAD NETSIM SPIKES
%
%	Loads spikes from a NETSIM binary file.
%
%	INPUT
%	filepath: /path/to/file
%
%	OUTPUT
%	ids - spike ids (int32) (note: 1 was added to this vector to adhere to Matlab's 1-based indexing)
%	times - spike times (double)
%

assert( ~isempty( regexp( filepath, 'spk\.bin?$', 'match', 'once' ) ) )

% open file
fid = fopen( filepath, 'rb' );

% read data
ids = fread( fid, inf, 'int32', 8 ); frewind( fid ); fseek(fid,4,'bof');
ids = ids + 1;
times = fread( fid, inf, 'double', 4 );

% close file
fclose( fid );
