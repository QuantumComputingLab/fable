function [UA, OA, alpha, info] = fable( A, compr_type, compr_val, logging )
% FABLE -- Fast Approximate BLock Encodings.
%
% INPUT
% -----
% A:            matrix to be block encoded
% compr_type:   type of compression algorithm:
%                 * 'percentage' : compr_val between 0-100%, x% largest
%                                  coefficients retained
%                 * 'cutoff'     : compr_vall determines cutoff value, larger
%                                  coefficients retained
% compr_val:    input parameter for compression algorithm
% logging:      true/false, if true info will log information about compression
%
% OUTPUT
% ------
% UA:       QCLAB circuit that block encodes A
% OA:       QCLAB matrix oracle circuit
% alpha:    subnormalization factor
% info:     struct containing some info on compression algorithm and circuit
%
% Copyright Daan Camps, Roel Van Beeumen, 2022.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % subnormalization
  alpha = norm(A(:), Inf) ;
  if alpha > 1
    A = A / alpha ;
  else
    alpha = 1 ;
  end
  N = size(A, 1) ;
  n = log2( N ) ;
  
  assert( N == 2^n ) ;
  assert( N == size(A,2) ) ;
  
  a = A(:) ;  
  
  % info struct
  if logging
    info = struct() ;
    info.vec.original = a ;
  else
    info = false ;
  end
  
  if isreal(A)
    if logging, info.datatype = 'real' ; end
    
    % Permuted FWHT
    a = 2 * acos( a ) ;
    a = sfwht( a ) ;
    a = grayPermutation( a ) ;
    
    if logging, info.vec.transformed = a ; end
    
    % Theshold vector according to compression criterion
    if strcmp( compr_type, 'percentage' )
      [~, sortIdx] = sort( abs(a),'ascend') ;
      cutoff = floor( (compr_val/100.0) * N^2 ) ;
      a( sortIdx(1:cutoff) ) = 0 ;
      if logging
        info.vec.compressed = a ;
        info.vec.zeroed = cutoff ;
      end
    elseif strcmp (compr_type, 'cutoff' )
      if logging, info.vec.zeroed = sum( abs(a) <= compr_val ) ; end
      a( abs(a) <= compr_val ) = 0 ;
      if logging, info.vec.compressed = a ; end
    end
    
    % Matrix oracle
    [OA, info_OA] = comprUniformRotation( a, n, N, 'RY', logging ) ;
    
    if logging
      info.circ.nCNOT = info_OA.nCNOT ;
      info.circ.nRY = info_OA.nG ;
    end
    
  else % complex case
    info.datatype = 'complex' ;
    
    % Generate magnitude and phase vectors
    a_m = zeros(size(a)) ;
    a_p = zeros(size(a)) ;
    for i = 1:length(a)
      if abs(imag(a(i))) < 1e-14
        % consider real (use sign in a_m )
        a_m(i) = real( a(i) ) ;
        a_p(i) = 0 ;
      else
        % consider complex ( a_m positive )
        a_m(i) = abs( a(i) ) ;
        a_p(i) = angle( a(i) ) ;
      end
    end
    
    if logging
      info.vec.magnitude.original = a_m ;
      info.vec.phase.original = a_p ;
    end
    
    % Permuted FWHT
    % magnitude
    a_m =  2 * acos( a_m ) ;
    a_m = sfwht( a_m ) ;
    a_m = grayPermutation( a_m ) ;
    % phase
    a_p = -2 * a_p ;
    a_p = sfwht( a_p ) ;
    a_p = grayPermutation( a_p ) ;
    
    if logging
      info.vec.magnitude.transformed = a_m ;
      info.vec.phase.transformed = a_p ;
    end
    
    % Theshold vectors according to compression criteria
    if strcmp( compr_type, 'percentage' )
      cutoff = floor( (compr_val/100.0) * N^2 ) ;
      [~, sortIdx] = sort( abs(a_m),'ascend') ;
      a_m( sortIdx(1:cutoff) ) = 0 ;
      [~, sortIdx] = sort( abs(a_p),'ascend') ;
      a_p( sortIdx(1:cutoff) ) = 0 ;
      if logging
        info.vec.magnitude.compressed = a_m ;
        info.vec.magnitude.zeroed = cutoff ;
        info.vec.phase.compressed = a_p ;
        info.vec.phase.zeroed = cutoff ;
      end
    elseif strcmp (compr_type, 'cutoff' )
      if logging 
        info.vec.magnitude.zeroed = sum( abs(a_m) <= compr_val ) ; 
        info.vec.phase.zeroed = sum( abs(a_p) <= compr_val ) ; 
      end
      a_m( abs(a_m) <= compr_val ) = 0 ;
      a_p( abs(a_p) <= compr_val ) = 0 ;
      if logging
        info.vec.magnitude.compressed = a_m ; 
        info.vec.phase.compressed = a_p ; 
      end
    end
    
    OA = qclab.QCircuit( 2*n + 1 ) ;
    [circ_mag, info_mag] = comprUniformRotation( a_m, n, N, 'RY', logging ) ;
    [circ_ph, info_ph] = comprUniformRotation( a_p, n, N, 'RZ', logging ) ;
    
    % Add both to OA   
    OA.push_back( circ_mag ) ;
    OA.push_back( circ_ph ) ;
    
    if logging
      info.circ.nCNOT = info_mag.nCNOT + info_ph.nCNOT;
      info.circ.nRY = info_mag.nG ;
      info.circ.nRZ = info_ph.nG ;
    end
  end

  % block-encoding circuit
  UA = qclab.QCircuit( 2*n + 1 ) ;
  
  % Diffusion on row indices
  for i = 1:n
    UA.push_back( qclab.qgates.Hadamard( i ) ) ;
  end
  
  % Matrix oracle
  UA.push_back( OA ) ;
  
  % SWAP register
  for i = 1:n
    UA.push_back( qclab.qgates.SWAP( i, 2*n - i + 1 ) ) ;
  end
  
  % Diffusion on row indices
  for i = 1:n
    UA.push_back( qclab.qgates.Hadamard( i ) ) ;
  end
  
  if logging
    info.circ.nH = 2*n ;
    info.circ.nSWAP = n ;
  end
    
end % end of fable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C, info] = comprUniformRotation( vec, n, N, type, logging )
% Compute a compressed uniform rotation circuit based on the thresholded vector
% vec of length N^2; N = 2^n, and type is either 'RY' or 'RZ'
  
  if strcmp( type, 'RY' )
    G = @qclab.qgates.RotationY ;
  elseif strcmp( type, 'RZ' )
    G = @qclab.qgates.RotationZ ;
  end
  
  C = qclab.QCircuit( 2*n + 1 ) ;
  
  nG = 0 ; nCNOT = 0 ;
  i = 1 ;
  while i <= N^2
    % Reset parity check
    parity_check = int32(0) ;
      
    % Add rotation gate
    if vec(i) ~= 0
      C.push_back( G( 0, vec(i) ) ) ;
      nG = nG + 1 ;
    end
      
    % Loop over sequence of consecutive zero angles
    while true
      % Compute control qubit
      ctrl = computeControl( i, n ) ;
      
      % update parity check
      if bitget( parity_check, ctrl ) 
        parity_check = bitset( parity_check, ctrl, int32(0) ) ;
      else
        parity_check = bitset( parity_check, ctrl, int32(1) ) ;
      end
        
      % update outer loop counter
      i = i + 1 ;
        
      % check exit condition
      if i > N^2 ||  vec(i) ~= 0
          break
      end
    end
      
    % Add CNOTs based on parity_check
    for j = 1:2*n
      if bitget( parity_check, j )
         C.push_back( qclab.qgates.CNOT( j, 0 ) ) ;
         nCNOT = nCNOT + 1 ;
      end
    end
      
  end
  
  if logging
    info = struct() ;
    info.nG = nG ;
    info.nCNOT = nCNOT ;
  else
    info = false ;
  end

end % end of comprUniformRotation

function [ctrl] = computeControl( i, n )
% Compute the control qubit based on the index i and the size n
  if i < 4^n
    ctrl = log2( bitxor( grayCode( i - 1 ), grayCode( i ) ) );
  else
    ctrl = 2*n - 1;
  end
  if ctrl > n - 1
    ctrl = 3*n - ctrl - 1 ;
  end
  ctrl = ctrl + 1 ;
end % end of computeControl

function x = grayCode( x )
  x = bitxor( x, bitshift( x, -1 ) );
end % end of grayCode

function [ b ] = grayPermutation( a )
  k = log2( size(a, 1) ) ; 
  b = zeros( 2^k, 1 );
  for i = 0 : 2^k - 1
    b( i + 1 ) = a( grayCode( i ) + 1 );
  end
end % end of grayPermutation

function [ a ] = sfwht( a )
% Scaled fast Walsh-Hadamard transform
  k = log2(size(a, 1) ) ;
  for h = 1:k
    for i = 1:2^h:2^k
      for j = i:i+2^(h-1)-1
        x = a( j );
        y = a( j + 2^(h-1) );
        a( j ) = ( x + y ) / 2 ;
        a( j + 2^(h-1) ) = ( x - y ) / 2 ;
      end
    end
  end
end % end of sfwht