function nbytes = sizeof(precision)
%SIZEOF  return the number of bytes of a builtin data type.
%   NBYTES = SIZEOF(PRECISION) returns the number of bytes of a single
%   element of class PRECISION.  PRECISION must be the name of one of the
%   builtin data types.
%
%   Knowing the number of bytes for a datatype is useful when performing
%   file I/O where some operations are defined in numbers of bytes.
%
%   Example:
%       nbytes = sizeof('single');


try
    z = zeros(1, precision); %#ok, we use 'z' by name later.
catch
    error('Unsupported class for finding size');
end

w = whos('z');
nbytes = w.bytes;

    
    
