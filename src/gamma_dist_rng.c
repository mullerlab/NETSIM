/********************************************************************************
 *                                                                              *
 *  FILE:   gamma_dist_rng.h                                                    *
 *  VERSION:  0.1                                                               *
 *  PROGRAM:  NETSIM                                                            *
 *                                                                              *
 *  PURPOSE:  Sampling from a gamma distribution    							*
 *																				*
 *  DETAILS:  Returns a deviate distributed as a gamma distribution of integer  *
 * 				order k, i.e., a waiting time to the kth event in a Poisson 	*	
 * 				process of unit mean, using UNIFORM_RAND as the source of 		*
 * 				uniform deviates. Adapted from Knuth, The Art of Computer  		*
 * 				Programming, Volume 2, 1997.									*
 *                                                                              *
 *  Copyright (C) 2016-2018 Lyle Muller                                         *
 *  http://cnl.salk.edu/~lmuller                                                *
 *                                                                              *
 * ---------------------------------------------------------------------------- *
 *                                                                              *
 *  DEVELOPMENT: Lyle Muller, Charlee Fletterman, Theo Desbordes                *
 *                                                                              *
 * ---------------------------------------------------------------------------- *
 *                                                                              *
 * This file is part of NETSIM.                                                 *
 *                                                                              *
 *     NETSIM is free software: you can redistribute it and/or modify           *
 *     it under the terms of the GNU General Public License as published by     *
 *     the Free Software Foundation, either version 3 of the License, or        *
 *     (at your option) any later version.                                      *
 *                                                                              *
 *     NETSIM is distributed in the hope that it will be useful,                *
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *     GNU General Public License for more details.                             *
 *                                                                              *
 *     You should have received a copy of the GNU General Public License        *
 *     along with NETSIM.  If not, see <http://www.gnu.org/licenses/>.          *
 *                                                                              *
 ********************************************************************************/

#include <stdio.h>
#include <math.h>

#include "rng.h"

/* function: gamma_sampling */
double gamma_sampling( int k, double scale, struct rng_state * rng )
{

	int ii;
	double x;

	if ( k < 1 ) { printf( "Error in routine gamma sampling: unsupported parameter k (shape of the gamma distribution)\n" ); }
	
	/* Use direct method, adding waiting times. */
	x = 1.0;
	for ( ii = 1 ; ii <= k ; ii++ ) x *= rng_dbl64(rng);
	x = -log(x);

	return x * scale;

}
