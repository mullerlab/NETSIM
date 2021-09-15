function [ data ] = load_netsim_data( filepath, varargin )
% *NETSIM*
%
%	LOAD NETSIM DATA
%
%	Loads data from a NETSIM binary file.
%
%	INPUT
%	filepath: /path/to/file
%   varargin: optional reshaping parameters
%
%	OUTPUT
%	data - data values (double)
%

assert( ~isempty( regexp( filepath, '\.bin?$', 'match', 'once' ) ) )

% open and read file
fid = fopen( filepath, 'rb' );
data = fread( fid, 'double' );

% optional reshape in function
if ( nargin > 1 ), data = reshape( data, varargin{1}, [] ); end
if ( nargin > 2 ), data = data(1:varargin{2},:); end

% close file
fclose( fid );
