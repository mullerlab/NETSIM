/********************************************************************************
 *                                                                              *
 *  FILE:   netsim.c                                                            *
 *  VERSION:  0.1                                                               *
 *  PROGRAM:  NETSIM                                                            *
 *                                                                              *
 *  PURPOSE:  A fast, large-scale simulator for topographic spiking networks    *
 *                                                                              *
 *  Copyright (C) 2016-2020 Lyle Muller                                         *
 *  http://mullerlab.ca                                                         *
 *                                                                              *
 * ---------------------------------------------------------------------------- *
 *                                                                              *
 *  DEVELOPMENT: Lyle Muller, Charlee Fletterman, Theo Desbordes, Gabriel       *
 *  Benigno, Christopher Steward                                                *
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


/* global includes */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <getopt.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <assert.h>

/* local includes */
#include "isaac64.h"
#include "gamma_dist_rng.h"
#include "rng.h"

/* defined constants */
#define NETSIM_NOERROR            0
#define NETSIM_ERROR              1
#define NETSIM_BAD_PARAMETER      2
#define NETSIM_EARLY_EXIT         -1

#define NETSIM_RANDOM_CONNECT     1
#define NETSIM_GAUSSIAN_CONNECT   2
#define NETSIM_TEST_CONNECT       3
#define NETSIM_GAUSSIAN_2D_CONNECT   4

#define NETSIM_INIT_STANDARD      0
#define NETSIM_INIT_SUSTAINED     1

#define FILE_BUFFER_SIZE          2048

/* individual parameter structures */
struct neuron_parameters {
  double taum, taur, tauref, vr, vth, vreset, Gl, El, Ie, Cm;
};

struct synapse_parameters {
  double ge, gi, taue, taui, Ee, Ei, synapse_delay, min_conduction_speed, max_conduction_speed, p_release, min_uniform_delay, max_uniform_delay;
};

struct network_parameters {
  int N, K;
  int spatial_dimensions;
  double L, sigma_space, poisson_rate, poisson_start_time, poisson_stop_time, poisEIratio, rewiring_probability;
};

struct simulation_parameters {

  /* simulation control parameters */
  int connector, initiator, output_file_code, bin_size, record_downsample_factor, buffer_length, cutoff_frequency;
  int conduction_speed_code, gamma_parameter_shape, output_path_override;
  double T, dt, vm_mean, vm_sigma, ge_mean, ge_sigma, gi_mean, gi_sigma, gamma_parameter_scale;
  char * parameter_file_path, * optional_connection_path, * output_path;

  /* recording parameters */
  int spikecount, save_connectivity;
  double start_record_time, stop_record_time, report_minutes;
  FILE * vfile, * spikefile, * gefile, * gifile;

  /* scan parameters */
  int job_id;

};

/* umbrella parameter structure */
struct parameters {
  struct neuron_parameters neuron;
  struct synapse_parameters synapse;
  struct network_parameters network;
  struct simulation_parameters simulation;
};

/* synapse bitfield structure */
struct trip {
  unsigned long offset : 64;
};

/* neuron structure */
struct neuron {
  bool excitatory;
  double v;
  double ge;
  double gi;
  struct trip * K;
  double conduction_speed;
  double x_position;
  double y_position;
  double lastSpike;
  int lastPoisE;
  int lastPoisI;
  double * buff_input_E;
  double * buff_input_I;
};

/* output structure */
struct output {
  int spikecount, spikecount_init;
  double endtime, vm_mean_all, p_release_out;
  double vmmean, vmsigma, gemean, gesigma, gimean, gisigma; /* end-state of the network */
  int NsideE; /* for 2d lattice (excitatory) */
  int NsideI; /* for 2d lattice (inhibitory) */
  FILE * vfile;
  FILE * spikefile;
  FILE * individual_traces;
  FILE * individual_traces_ge;
  FILE * gefile;
  FILE * gifile;
};

#include "helper_functions/helper_scan.h"
#include "helper_functions/helper_functions.h"
#include "helper_functions/helper_gaussian.h"

/* gaussian connect */
int gaussian_connect( struct network_parameters * net, struct simulation_parameters * sim, struct synapse_parameters * syn, struct neuron * network, struct rng_state * rng )
{

  /* init */
  /* special for calculating netsim $% */
  bool *temparr;

  /* bitfield variables */
  int trip_number, trip_offset;
  const unsigned long mask[3] = {0x7FFFFFFFFFEUL, 0xFFFFF800003FFFFEUL, 0xFFFFFFFFFFC00000UL};
  const unsigned short bitshift[3] = {43,22,1};
  struct trip * tmp;

  /* same for calculating and saved netsim */
  int ii, jj;
  int Ne, Ni;
  int count, count_connections;
  double x_position, distance, dx_exc, dx_inh;

  double max_distance = 0;

  /* calculate */
  Ne = round( .8*net->N );
  Ni = round( .2*net->N );
  dx_exc = net->L / Ne;
  dx_inh = net->L / Ni;

  /* temporarily save connections in boolean size N, then translate to sparse */
  /* like this it is way faster to check for double connections */
  temparr = (bool *) malloc( net->N * sizeof(bool) ); /* $% */

  for ( ii = 0 ; ii < net->N ; ii++ )
  {

    x_position = (*(network+ii)).x_position;

    /* reset temparr */
    memset( temparr, (bool) 0, (size_t) net->N * sizeof(bool) );

    count = (int) round(.8*net->K);
    count_connections = 0;
    while ( count > 0 )
    {

      /* draw gaussian distance */
      distance = rng_gauss( rng ) * net->sigma_space;

      /* within range? */
      if ( fabs(distance) > (.5 * net->L) ) {  continue;  }

      jj = calculate_closest( distance, x_position, dx_exc , Ne );

      /* no self- or double connections */
      if ( ((ii==jj) && (*(network+ii)).excitatory) || *(temparr+jj) ) {  continue;  }

      /* check for the largest connection distance */
      if ( fabs(distance) > max_distance ) { max_distance = fabs(distance); }

      trip_number = (int) floor( (double) count_connections / (double) 3 );
      trip_offset = ( count_connections % 3 );

      tmp = ( ( ( *(network+ii) ).K ) + trip_number );
      tmp->offset = ( ( tmp->offset & mask[trip_offset] ) ^ ( ( (unsigned long) jj ) << bitshift[trip_offset] ) );

      *(temparr+jj) = 1; /* $% */
      count--;
      count_connections++;

    }

    count = (int) round(.2*net->K);
    while ( count > 0 )
    {

      distance = rng_gauss( rng) * net->sigma_space;

      if ( fabs(distance) > (.5 * net->L) ) {  continue;  }

      jj = Ne + calculate_closest( distance, x_position, dx_inh, Ni );

      /* no self- or double connections */
      if ( ((ii==jj) && !(*(network+ii)).excitatory) || *(temparr+jj) ) {  continue;  }

      /* check for the largest connection distance */
      if ( fabs(distance) > max_distance ) { max_distance = fabs(distance); }
      
      trip_number = (int) floor( (double) count_connections / (double) 3 );
      trip_offset = ( count_connections % 3 );

      tmp = ( ( ( *(network+ii) ).K ) + trip_number );
      tmp->offset = ( ( tmp->offset & mask[trip_offset] ) ^ ( ( (unsigned long) jj ) << bitshift[trip_offset] ) );

      *(temparr+jj) = 1; /* $% */
      count--;
      count_connections++;

    }

  }

  /* assign buffer length based on max distance */
  sim->buffer_length = (int) round( ( syn->synapse_delay + (max_distance / (syn->min_conduction_speed)) ) / sim->dt ) + 1; 

  /* cleanup */
  free( temparr );
  return NETSIM_NOERROR;

}

/* gaussian connect */
int gaussian_connect_2D( struct network_parameters * net, struct simulation_parameters * sim, struct synapse_parameters * syn, struct neuron * network, struct rng_state * rng )
{

  /* init */
  /* special for calculating netsim $% */
  bool *temparr;

  /* bitfield variables */
  int trip_number, trip_offset;
  const unsigned long mask[3] = {0x7FFFFFFFFFEUL, 0xFFFFF800003FFFFEUL, 0xFFFFFFFFFFC00000UL};
  const unsigned short bitshift[3] = {43,22,1};
  struct trip * tmp;

  /* same for calculating and saved netsim */
  int ii, jj, c1, c2;
  int Ne, Ni, NrowE, NrowI;
  int count, count_connections;
  double x_position, y_position, distance_x, distance_y, distance, dx_exc, dx_inh;

  double max_distance = 0;

  /* calculate */
  Ne = round( .8*net->N ); Ni = round( .2*net->N );
  NrowE = floor( sqrt(Ne) ); NrowI = floor( sqrt(Ni) );
  dx_exc = net->L / NrowE; dx_inh = net->L / NrowI;

  /* temporarily save connections in boolean size N, then translate to sparse */
  /* like this it is way faster to check for double connections */
  temparr = (bool *) malloc( net->N * sizeof(bool) ); /* $% */

  for ( ii = 0 ; ii < net->N ; ii++ )
  {

    x_position = (*(network+ii)).x_position; y_position = (*(network+ii)).y_position;

    /* reset temparr */
    memset( temparr, (bool) 0, (size_t) net->N * sizeof(bool) );

    count = (int) round(.8*net->K);
    count_connections = 0;
    while ( count > 0 )
    {

      /* draw gaussian distance */
      distance_x = rng_gauss( rng ) * net->sigma_space;
      distance_y = rng_gauss( rng ) * net->sigma_space;
      distance = sqrt( pow(distance_x,2) + pow(distance_y,2) );

      /* within range? */
      if ( fabs(distance) > (M_SQRT2 * .5 * net->L) ) {  continue;  }

      c1 = calculate_closest( distance_x, x_position, dx_exc, NrowE ); /* x-coordinate */
      c2 = calculate_closest( distance_y, y_position, dx_exc, NrowE ); /* y-coordinate */
      jj = c1 + ( c2 * NrowE ); 

      /* no self- or double connections */
      if ( ((ii==jj) && (*(network+ii)).excitatory) || *(temparr+jj) ) {  continue;  }

      /* check for the largest connection distance */
      if ( fabs(distance) > max_distance ) {  max_distance = fabs(distance); }

      trip_number = (int) floor( (double) count_connections / (double) 3 );
      trip_offset = ( count_connections % 3 );

      tmp = ( ( ( *(network+ii) ).K ) + trip_number );
      tmp->offset = ( ( tmp->offset & mask[trip_offset] ) ^ ( ( (unsigned long) jj ) << bitshift[trip_offset] ) );

      *(temparr+jj) = 1; /* $% */
      count--; count_connections++;

    }

    count = (int) round(.2*net->K);
    while ( count > 0  )
    {

      /* draw gaussian distance */
      distance_x = rng_gauss( rng ) * net->sigma_space;
      distance_y = rng_gauss( rng ) * net->sigma_space;
      distance = sqrt( pow(distance_x,2) + pow(distance_y,2) );

      if ( fabs(distance) > (M_SQRT2 * .5 * net->L) ) {  continue;  }

      c1 = calculate_closest( distance_x, x_position, dx_inh, NrowI ); /* x-coordinate */
      c2 = calculate_closest( distance_y, y_position, dx_inh, NrowI ); /* y-coordinate */
      jj = Ne + c1 + ( c2 * NrowI );

      /* no self- or double connections */
      if ( ((ii==jj) && !(*(network+ii)).excitatory) || *(temparr+jj) ) {  continue;  }

      /* check for the largest connection distance */
      if ( fabs(distance) > max_distance ) { max_distance = fabs(distance); }
      
      trip_number = (int) floor( (double) count_connections / (double) 3 );
      trip_offset = ( count_connections % 3 );

      tmp = ( ( ( *(network+ii) ).K ) + trip_number );
      tmp->offset = ( ( tmp->offset & mask[trip_offset] ) ^ ( ( (unsigned long) jj ) << bitshift[trip_offset] ) );

      *(temparr+jj) = 1; /* $% */
      count--;
      count_connections++;

    }

  }

  /* assign buffer length based on max distance */
  sim->buffer_length = (int) round( ( syn->synapse_delay + (max_distance / (syn->min_conduction_speed)) ) / sim->dt ) + 1; 

  /* cleanup */
  free( temparr );
  return NETSIM_NOERROR;

}

/* function: Erdos-Renyi random graph */
int random_connect( struct network_parameters * net, struct neuron * network, struct rng_state * rng )
{

  /* init */
  register int ii, jj;

  /* bitfield variables */
  const unsigned long mask[3] = {0x7FFFFFFFFFEUL, 0xFFFFF800003FFFFEUL, 0xFFFFFFFFFFC00000UL};
  const unsigned short bitshift[3] = {43,22,1};
  uint32_t target = 0;
  struct trip * tmp;

  /* loop over network */
  for ( ii = 0 ; ii < net->N ; ii++ )
  {

    for ( jj = 0 ; jj < floor( (double) net->K / (double) 3 ) ; jj++ )
    {

      tmp = ( ( ( *(network+ii) ).K ) + jj );

      target = rng_uint32( rng ) % net->N;
      tmp->offset = ( ( tmp->offset & mask[0] ) ^ ( ( (unsigned long) target ) << bitshift[0] ) );

      target = rng_uint32( rng ) % net->N;
      tmp->offset = ( ( tmp->offset & mask[1] ) ^ ( ( (unsigned long) target ) << bitshift[1] ) );

      target = rng_uint32( rng ) % net->N;
      tmp->offset = ( ( tmp->offset & mask[2] ) ^ ( ( (unsigned long) target ) << bitshift[2] ) );

    }

    for ( jj = 0 ; jj < ( net->K % 3 ) ; jj++ )
    {

      tmp = ( ( ( *(network+ii) ).K ) + (int) floor( (double) net->K / (double) 3 ) );

      target = rng_uint32( rng ) % net->N;
      tmp->offset = ( ( tmp->offset & mask[jj] ) ^ ( ( (unsigned long) target ) << bitshift[jj] ) );      

    }

  }

  return NETSIM_NOERROR;

}

/* function: load test conections from file */
int test_connect( struct parameters * p, struct neuron * network )
{

  /* init */
  int ii, jj, id;
  FILE * file;
  file = fopen( p->simulation.optional_connection_path, "r" );

  int trip_number, trip_offset;
  const unsigned long mask[3] = {0x7FFFFFFFFFEUL, 0xFFFFF800003FFFFEUL, 0xFFFFFFFFFFC00000UL};
  const unsigned short bitshift[3] = {43,22,1};
  struct trip * tmp;
  
  /* loop over network */
  for ( ii = 0 ; ii < p->network.N ; ii++ )
  {

    for ( jj = 0 ; jj < p->network.K ; jj++ ) 
    {  
      
      /* get connection from file */
      fscanf( file, "%d\n", &id );  

      /* write connection to neuron structure */
      trip_number = (int) floor( (double) jj / (double) 3 );
      trip_offset = ( jj % 3 );

      tmp = ( ( ( *(network+ii) ).K ) + trip_number );
      tmp->offset = ( ( tmp->offset & mask[trip_offset] ) ^ ( ( (unsigned long) id ) << bitshift[trip_offset] ) );

    }

  }

  /* clean up */
  fclose( file );
  return 0;

}

/* function: random rewiring of EXCITATORY network connections */
int random_rewiring( struct parameters * p, struct neuron * network, struct rng_state * rng )
{

  /* init */
  int ii, jj, Ne;

  int trip_number, trip_offset;
  const unsigned long mask[3] = {0x7FFFFFFFFFEUL, 0xFFFFF800003FFFFEUL, 0xFFFFFFFFFFC00000UL};
  const unsigned short bitshift[3] = {43,22,1};
  uint32_t target = 0;
  struct trip * tmp;
  
  /* number of excitatory neurons */
  Ne = round( .8*p->network.N );

  /* loop over network */
  for ( ii = 0 ; ii < Ne ; ii++ )
  {

    for ( jj = 0 ; jj < p->network.K ; jj++ ) 
    {

      if ( rng_dbl64( rng ) < p->network.rewiring_probability )  
      {

        /* write connection to neuron structure */
        trip_number = (int) floor( (double) jj / (double) 3 );
        trip_offset = ( jj % 3 );

        /* select random target */
        target = rng_uint32( rng ) % p->network.N;

        /* write random target */
        tmp = ( ( ( *(network+ii) ).K ) + trip_number );
        tmp->offset = ( ( tmp->offset & mask[trip_offset] ) ^ ( ( (unsigned long) target ) << bitshift[trip_offset] ) );

      }
      
    }

  }

  /* assign buffer length as in random graph simulation */
  p->simulation.buffer_length = (int) ( ( p->synapse.synapse_delay + ((.5 * p->network.L) / (p->synapse.min_conduction_speed)) ) / p->simulation.dt ) + 1;     

  /* clean up */
  return 0;

}

/* function: init to zero */
int initialize_simulation( struct parameters * p, struct neuron * network, struct output * o, struct rng_state * rng )
{

  /* init */
  register int ii;
  int Ne = (int) round( .8 * p->network.N ) ;
  int buffer_length;
  double dx_exc, dx_inh;
  double tmp_conduction_speed;
  char str[250];

  /* calculate */
  dx_exc = p->network.L / Ne;
  dx_inh = p->network.L / ( .2 * p->network.N );

  /* calculate needed buffer size */
  if ( p->synapse.max_uniform_delay > 0 ) {  buffer_length = (int) ( p->synapse.max_uniform_delay / p->simulation.dt ) + 1;  }
  else {  buffer_length = (int) ( ( p->synapse.synapse_delay + ((.5 * p->network.L) / (p->synapse.min_conduction_speed)) ) / p->simulation.dt ) + 1;  }
  p->simulation.buffer_length = buffer_length;

  /* initialize network */
  for ( ii = 0 ; ii < p->network.N ; ii++ )
  {

    ( *(network+ii) ).excitatory    = ( ii < Ne );
    ( *(network+ii) ).v             = p->neuron.El;
    ( *(network+ii) ).ge            = 0;
    ( *(network+ii) ).gi            = 0;
    ( *(network+ii) ).K             = (struct trip *) malloc( sizeof(struct trip) * ceil( (double) p->network.K / (double) 3 ) );
    ( *(network+ii) ).lastSpike     = -1;
    ( *(network+ii) ).lastPoisE     = -1;
    ( *(network+ii) ).lastPoisI     = -1;

    /* calculate neuron's x_position in the ring */
    if ( ( *(network+ii) ).excitatory == 1 ) {  (*(network+ii) ).x_position = ii * dx_exc;  }
    else {  ( *(network+ii) ).x_position = ( ii - Ne ) * dx_inh;  }
    ( *(network+ii) ).y_position = 0;

    /* calculate conduction velocity: code 0 -> constant speed, 1 -> uniform, 2 -> gamma distribution  */
    if ( p->simulation.conduction_speed_code == 0 ) {  ( *(network+ii) ).conduction_speed = p->synapse.min_conduction_speed;  }
    else if ( p->simulation.conduction_speed_code == 1 ) {  ( *(network+ii) ).conduction_speed =
      ( p->synapse.max_conduction_speed - p->synapse.min_conduction_speed ) * rng_dbl64(rng) + p->synapse.min_conduction_speed;  }
    else if ( p->simulation.conduction_speed_code == 2 ) {  

      do { tmp_conduction_speed = gamma_sampling( p->simulation.gamma_parameter_shape, p->simulation.gamma_parameter_scale, rng ); }
      while ( (tmp_conduction_speed < p->synapse.min_conduction_speed) || (tmp_conduction_speed > p->synapse.max_conduction_speed) );

      ( *(network+ii) ).conduction_speed = tmp_conduction_speed;

    }

    assert( ( *(network+ii) ).K!=NULL );

  }

  /* initialize output files */
  /* REMEMBER WHEN ADDING FILES, TO MOVE THEM FROM SCRATCH IF RUN WITH JOBMANAGER! */
  if ( p->simulation.stop_record_time > 0 )
  {

    sprintf( str,  "%s/%08dvms.bin", p->simulation.output_path, p->simulation.output_file_code );
    o->vfile = fopen( str, "w" );
    sprintf( str, "%s/%08dge.bin", p->simulation.output_path, p->simulation.output_file_code );
    o->gefile = fopen( str, "w" );
    sprintf( str, "%s/%08dgi.bin", p->simulation.output_path, p->simulation.output_file_code );
    o->gifile = fopen( str, "w" );

  }
  /* if stop_record_time == 0 nothing will be recorded so no point in making the file */
  else
  {

    o->vfile = fopen( "/dev/null", "w" );
    o->gefile = fopen( "/dev/null", "w" );
    o->gifile = fopen( "/dev/null", "w" );

  }

  sprintf( str,  "%s/%08dspk.bin", p->simulation.output_path, p->simulation.output_file_code );
  o->spikefile = fopen( str, "wb" );
  sprintf( str,  "%s/%08dindividual.bin", p->simulation.output_path, p->simulation.output_file_code );
  o->individual_traces = fopen( str, "wb" );
  sprintf( str,  "%s/%08dindividualge.bin", p->simulation.output_path, p->simulation.output_file_code );
  o->individual_traces_ge = fopen( str, "wb" );

  assert( o->vfile!=NULL && o->spikefile!=NULL && o->individual_traces!=NULL );

  return 0;

}

/* function: allocate input buffers */
int allocate_input_buffers( struct parameters * p, struct neuron * network )
{

  /* init */
  int ii, buffer_length; 
  buffer_length = p->simulation.buffer_length;

  /* initialize buffers in network */
  for ( ii = 0 ; ii < p->network.N ; ii++ )
  {

    ( *(network+ii) ).buff_input_E  = (double *) calloc( buffer_length, sizeof(double) );
    ( *(network+ii) ).buff_input_I  = (double *) calloc( buffer_length, sizeof(double) );

    /* check memory allocation */
    assert( ( *(network+ii) ).buff_input_E!=NULL && ( *(network+ii) ).buff_input_I!=NULL );

  }

  return NETSIM_NOERROR;

}

/* function: intialize for self-sustained activity */
int initialize_sustained_simulation( struct parameters * p, struct neuron * network, struct output * o, struct rng_state * rng )
{

  printf( "> Initializing self-sustained activity " );

  if ( p->simulation.initiator == 0 )
  {

    printf( "with gaussian vm, ge, and gi\n" );

    /* init */

    register int ii;

    for ( ii = 0 ; ii < p->network.N ; ii++ )
    {
      /* init v, ge and gi with normal distribution for sustained (emperical) */
      ( *(network+ii) ).v             = rng_gauss( rng ) * p->simulation.vm_sigma + p->simulation.vm_mean;
      ( *(network+ii) ).ge            = rng_gauss( rng ) * p->simulation.ge_sigma + p->simulation.ge_mean;
      ( *(network+ii) ).gi            = rng_gauss( rng ) * p->simulation.gi_sigma + p->simulation.gi_mean;
    }

  }
  else if ( p->simulation.initiator == 1 )
  {

    printf( "with poisson input\n" );

    /* init */
    int ii, jj;
    double hold_poisEIratio, hold_poisson_rate;
    double hold_poisson_start_time, hold_poisson_stop_time, hold_T;
    FILE * hold_spikefile;
    FILE * hold_vfile;
    FILE * hold_gefile;
    FILE * hold_gifile;
    FILE * hold_individual_traces;
    FILE * hold_individual_traces_ge;
    int hold_spikecount;
    int buffer_length, buffer_position;
    char str[2048];

    /* calculate needed buffer size */
    buffer_length = p->simulation.buffer_length;
    double * temparrE = calloc( buffer_length, sizeof(double) );
    double * temparrI = calloc( buffer_length, sizeof(double) );

    /* save for later to conserve original settings */
    hold_poisson_rate = p->network.poisson_rate;
    hold_poisEIratio = p->network.poisEIratio;
    hold_poisson_start_time = p->network.poisson_start_time;
    hold_poisson_stop_time = p->network.poisson_stop_time;
    hold_spikefile = o->spikefile;
    hold_vfile = o->vfile;
    hold_gefile = o->gefile;
    hold_gifile = o->gifile;
    hold_individual_traces = o->individual_traces;
    hold_individual_traces_ge = o->individual_traces_ge;
    hold_spikecount = o->spikecount;
    hold_T = p->simulation.T;

    /* set */
    p->network.poisson_rate = 2000 / p->network.K;
    p->network.poisEIratio = 1;
    p->network.poisson_start_time = 0;
    p->network.poisson_stop_time = 200e-3;
    sprintf( str,  "%s/%08dspk_init.bin", p->simulation.output_path, p->simulation.output_file_code );
    o->spikefile = fopen( str, "wb" ); /* reset this file every time initializer runs */
    p->simulation.T = 300e-3;
    o->vfile = fopen( "/dev/null", "rb" );
    o->gefile = fopen( "/dev/null", "rb" );
    o->gifile = fopen( "/dev/null", "rb" );
    o->individual_traces = fopen( "/dev/null", "rb" );
    o->individual_traces_ge = fopen( "/dev/null", "rb" );

    /* consistency check */
    if ( p->synapse.p_release > 0 && p->synapse.p_release < 1)
    {  printf( "Currently no external_input with p_release implemented. Exiting..\n" ); return -1;  }

    fclose( o->spikefile );

    buffer_position = (int) round( o->endtime / p->simulation.dt ) % buffer_length;

    for ( ii=0; ii<p->network.N; ii++ )
    {

      (*(network+ii)).lastSpike = (*(network+ii)).lastSpike - p->simulation.T;
      (*(network+ii)).lastPoisE = (*(network+ii)).lastPoisE - (int) round( p->simulation.T/p->simulation.dt );
      (*(network+ii)).lastPoisI = (*(network+ii)).lastPoisI - (int) round( p->simulation.T/p->simulation.dt );

      for ( jj=0; jj<buffer_length; jj++ )
      {
        *(temparrE+jj) = *( (*(network+ii)).buff_input_E + ( buffer_position + jj ) % buffer_length );
        *(temparrI+jj) = *( (*(network+ii)).buff_input_I + ( buffer_position + jj ) % buffer_length );
      }

      for ( jj=0; jj<buffer_length; jj++ )
      {
        *( (*(network+ii)).buff_input_E + jj) = *(temparrE+jj);
        *( (*(network+ii)).buff_input_I + jj) = *(temparrI+jj);
      }

    }

    p->network.poisson_rate = hold_poisson_rate;
    p->network.poisEIratio = hold_poisEIratio;
    p->network.poisson_start_time = hold_poisson_start_time;
    p->network.poisson_stop_time = hold_poisson_stop_time;
    p->simulation.T = hold_T;
    o->spikefile = hold_spikefile;
    o->spikecount = hold_spikecount;
    o->vfile = hold_vfile;
    o->gefile = hold_gefile;
    o->gifile = hold_gifile;
    o->individual_traces = hold_individual_traces;
    o->individual_traces_ge = hold_individual_traces_ge;

  }

  return 0;

}

/* function initialize 2D simulation */
int initialize_2D_simulation( struct parameters * p, struct neuron * network )
{

  /* init */
  int Ne = (int) round( .8 * p->network.N ); 
  int Ni = p->network.N - Ne;
  register int ii, jj;
  int NrowE, NrowI, index; 

  double dx_exc, dx_inh;

  /* calculate */
  NrowE = floor( sqrt(Ne) ); NrowI = floor( sqrt(Ni) );
  dx_exc = p->network.L / NrowE; dx_inh = p->network.L / NrowI;

  /* loop over neurons and calculate lattice indices */
  for ( ii = 0 ; ii < NrowE ; ii++ )
  {
    for ( jj = 0 ; jj < NrowE ; jj++ )
    {
      /* calculate index */
      index = jj + ( ii * NrowE );

      ( *(network+index) ).x_position = jj * dx_exc; 
      ( *(network+index) ).y_position = ii * dx_exc;
    }
  }

  /* loop over neurons and calculate lattice indices */
  for ( ii = 0 ; ii < NrowI ; ii++ )
  {
    for ( jj = 0 ; jj < NrowI ; jj++ )
    {
      /* calculate index */
      index = jj + ( ii * NrowI ) + Ne;

      ( *(network+index) ).x_position = jj * dx_inh; 
      ( *(network+index) ).y_position = ii * dx_inh;
    }
  }

  return NETSIM_NOERROR;

}

/* function: run_simulation */
int run_simulation( struct parameters * p, struct neuron * network, struct output * o, struct rng_state * rng )
{

  /* make local structure copies */
  struct network_parameters net = p->network;
  struct neuron_parameters neu = p->neuron;
  struct synapse_parameters syn = p->synapse;
  struct simulation_parameters sim = p->simulation;

  /* init local variables */
  register int ii, jj, kk;
  int tt, N, K, Nsteps;
  int buffer_length, buffer_position, insert_position;
  int number_of_seconds_between_progress_report;
  int spikecount_first_timecheck;
  double t, dt, time_first_timecheck, rate;
  double taue, taui, taur, taum, Ee, Ei, vreset, vth, ge, gi, Gl, El, Ie, Cm;
  clock_t begin, intermediate;

  /* bitfield variables */
  const unsigned long mask[3] = {0xFFFFF80000000000UL, 0x7FFFFC00000UL, 0x3FFFFEUL};
  const unsigned short bitshift[3] = {43,22,1};
  uint32_t target = 0;
  struct trip * tmp;

  /* init special for recording ^& */
  double * record_1_timestep, * record_1_timestep_binned;
  double * record_1_timestep_ge, * record_1_timestep_binned_ge;
  double * record_1_timestep_gi, * record_1_timestep_binned_gi;
  double totalsum_v, sum_v, sum_ge, sum_gi;
  char * spike_array;
  char * start_spike_array;
  int spike_index, vsteps, nvsteps;
  uint32_t origin;
  int number_of_recorded_individuals, record_individuals_jump;
  int bin_size, record_downsample_factor, record_start_bin, record_stop_bin;
  int spikes_in_this_round;
  unsigned long long totalcount_v;
  int NrowE;

  /* special init for p relase */
  double p_release;
  unsigned int spike_transmitted, spike_not_transmitted;
  spike_transmitted = spike_not_transmitted = 0;

  /* special init for recalculating delay */
  int delay_in_bins;
  double distance_x, distance_y, distance, L, synapse_delay;

  /* make local variable copies */
  N = net.N; K = net.K; L = net.L;
  taur = neu.taur; taum = neu.taum; vreset = neu.vreset; vth = neu.vth; El = neu.El; Ie = neu.Ie;
  taue = syn.taue; taui = syn.taui; Ee = syn.Ee; Ei = syn.Ei; ge = syn.ge; gi = syn.gi; synapse_delay = syn.synapse_delay; p_release = syn.p_release;
  dt = sim.dt; bin_size = sim.bin_size; record_downsample_factor = sim.record_downsample_factor;

  #ifdef EXTERNAL_INPUT
  /* start POISSON */
  int N_poisson_E, N_poisson_I, delay_in_bins_poisson;
  int stop_poisson_bin, start_poisson_bin;
  int count_poisE, count_poisI;
  double poisson_rate, poisson_ISI;
  poisson_rate = net.poisson_rate;
  start_poisson_bin = (int) round( net.poisson_start_time / sim.dt );
  stop_poisson_bin = (int) round( net.poisson_stop_time / sim.dt );
  N_poisson_E = (int) round( net.poisEIratio * net.K );
  N_poisson_I = net.K - N_poisson_E;
  count_poisE = 0; count_poisI = 0;
  /* end POISSON */
  #endif

  /* init recording */
  number_of_recorded_individuals = 10;
  record_individuals_jump = N / number_of_recorded_individuals;

  if ( p->network.spatial_dimensions == 1 ) {  vsteps = (int) floor( N / bin_size ); }
  if ( p->network.spatial_dimensions == 2 ) 
  	{  vsteps = (int) floor( floor(0.8*N) / (bin_size*bin_size) ); nvsteps = (int) floor( sqrt(vsteps) );  }
  NrowE = floor( sqrt(floor(0.8*N)) );

  record_start_bin = (int) floor( sim.start_record_time / dt );
  record_stop_bin = (int) floor( sim.stop_record_time / dt );
  record_1_timestep = (double *) malloc ( sizeof(double) * N );
  record_1_timestep_binned = (double *) malloc ( sizeof(double) * (int) vsteps );
  record_1_timestep_ge = (double *) malloc ( sizeof(double) * N );
  record_1_timestep_binned_ge = (double *) malloc ( sizeof(double) * (int) vsteps );
  record_1_timestep_gi = (double *) malloc ( sizeof(double) * N );
  record_1_timestep_binned_gi = (double *) malloc ( sizeof(double) * (int) vsteps );
  spike_array = (char *)malloc( ( sizeof(int) + sizeof(double) ) * N );

  start_spike_array = spike_array;
  fwrite( &number_of_recorded_individuals, sizeof(int), 1, o->individual_traces );

  /* calculate parameters */
  Cm = neu.Cm; Gl = Cm / taum;
  Nsteps = round( sim.T / dt );
  number_of_seconds_between_progress_report = sim.report_minutes * 60;
  buffer_length = sim.buffer_length;

  /* initialize variables and arrays */
  begin = clock(); intermediate = clock();
  spikes_in_this_round = 0; totalsum_v = 0; totalcount_v = 0;
  spikecount_first_timecheck = 0; time_first_timecheck = 0;

  /* MAIN SIMULATION LOOP */
  for ( tt = 0 ; tt < Nsteps ; tt++ )
  {

    t = dt * tt;
    buffer_position = tt % buffer_length;

    spike_array = start_spike_array; /* ^& */
    spikes_in_this_round = 0; /* ^& */

    if ( ( time_first_timecheck == 0 ) && ( t > 5e-3 ) )
    {
      spikecount_first_timecheck = o->spikecount; time_first_timecheck = t;
    }

    /* near to no slow-down (between 4e-5 and 5e-4 sec per timestep) */
    if ( ( ((int) ( clock() - intermediate ) / CLOCKS_PER_SEC ) > number_of_seconds_between_progress_report ) && ( t > 10e-3 ) )
    {

      rate = (double) ( o->spikecount - spikecount_first_timecheck ) / N / ( t - time_first_timecheck );
      printf("%f out of %f seconds simulated in %f (number of spikes so far: %d -> %f (Hz) )\n",
        t, sim.T, (double)( clock() - begin)/ CLOCKS_PER_SEC, o->spikecount, rate );
      fflush( stdout );
      intermediate = clock();

      /* if the rate is too high, stop the simulation */
      if ( rate > sim.cutoff_frequency )
      {
          printf( "Early exit because firing rate > %d Hz \n", sim.cutoff_frequency);
          o->endtime = t;
          o->vm_mean_all = totalsum_v / totalcount_v;

          return NETSIM_EARLY_EXIT;
      }

    }

    for ( ii = 0 ; ii < N ; ii++ )
    {

      #ifdef EXTERNAL_INPUT
      /* start POISSON */
      if ( ( (tt >= start_poisson_bin) && (tt < stop_poisson_bin) ) )
      {

        if ( ( N_poisson_E==0) || (*(network+ii)).lastPoisE > tt ) {  ;  }
        else
        {

          while ( (*(network+ii)).lastPoisE <= tt )
          {

            if ( (*(network+ii)).lastPoisE == tt ) /* need to do this for edge case of start, not too happy with it */
            {
              *(( *(network + ii)).buff_input_E + buffer_position ) += ge;
              count_poisE++;
            }

            poisson_ISI = ( -log( rng_dbl64( rng ) ) / ( poisson_rate * N_poisson_E ) );
            delay_in_bins_poisson = (int) round( poisson_ISI / dt );
            ( *(network+ii) ).lastPoisE = tt + delay_in_bins_poisson;

          }

        }

        if ( ( N_poisson_I==0) || (*(network+ii)).lastPoisI > tt ) {  ;  }
        else
        {

          while ( ( (*(network+ii)).lastPoisI <= tt  ) )
          {

            if ( (*(network+ii)).lastPoisI == tt ) /* need to do this for edge case of start, not too happy with it */
            {
              *(( *(network + ii)).buff_input_I + buffer_position ) += gi;
              count_poisI++;
            }

            poisson_ISI =  ( -log( rng_dbl64( rng ) ) / ( poisson_rate * N_poisson_I ) );
            delay_in_bins_poisson = (int) round( poisson_ISI / dt );
            ( *(network+ii) ).lastPoisI = tt + delay_in_bins_poisson;

          }

        }

      }
      /* end POISSON */
      #endif

      /* integrate synaptic input */

      ( *(network+ii) ). ge +=  -dt*( ( *(network+ii) ).ge )/taue +
                                *(( *(network+ii) ).buff_input_E + buffer_position);
      ( *(network+ii) ). gi +=  -dt*( ( *(network+ii) ).gi )/taui +
                                *(( *(network+ii) ).buff_input_I + buffer_position);


      /* reset buffer */
      *( ( *(network+ii) ).buff_input_E + buffer_position ) = 0;
      *( ( *(network+ii) ).buff_input_I + buffer_position ) = 0;

      if ( t > ( (*(network+ii)).lastSpike + taur ) )
      {

        ( *(network+ii) ).v +=
                dt * (
                  ( (*(network+ii)).ge*(Ee-(*(network+ii)).v) +  /* excitatory synapse */
                    (*(network+ii)).gi*(Ei-(*(network+ii)).v) +  /* inhibitory synapse */
                    Gl * (El - (*(network+ii)).v ) +             /* leak conductance */
                    Ie ) / Cm ) ;                                /* current injection */

        /* save membrane potential into temporary array ^& */
        *( record_1_timestep + ii ) = (*(network+ii)).v;
        *( record_1_timestep_ge + ii ) = (*(network+ii)).ge;
        *( record_1_timestep_gi + ii ) = (*(network+ii)).gi;

        if ( ( *(network+ii) ).v >= vth )
        {

          /* increment spike counter ^& */
          o->spikecount++;

          /* complicated way to save the spike index and spike time in a char array ^& */
          spike_index = ii;
          memcpy( (void *)spike_array, (void *)(&spike_index), sizeof(int) ); spike_array += sizeof(int);
          memcpy( (void *)spike_array, (void *)(&t), sizeof(double) ); spike_array += sizeof(double);
          spikes_in_this_round++;
          /* end ^& */

          ( *(network+ii) ).v = vreset;
          ( *(network+ii) ).lastSpike = t;

          if ( ( *(network+ii) ).excitatory )
          {

            for ( jj = 0 ; jj < floor( (double) K / (double) 3 ) ; jj ++ )
            {

              tmp = ( (*(network+ii)).K + jj );

              #include "templates/release_probability_template_start.h"

              target = (uint32_t) ( ( tmp->offset & mask[0] ) >> bitshift[0] );
              #include "templates/distance_template.h"
              *( (*(network + target)).buff_input_E + insert_position ) += ge;

              #include "templates/release_probability_template_stop.h"

              #include "templates/release_probability_template_start.h"

              target = (uint32_t) ( ( tmp->offset & mask[1] ) >> bitshift[1] );
              #include "templates/distance_template.h"
              *( (*(network + target)).buff_input_E + insert_position ) += ge;

              #include "templates/release_probability_template_stop.h"

              #include "templates/release_probability_template_start.h"

              target = (uint32_t) ( ( tmp->offset & mask[2] ) >> bitshift[2] );
              #include "templates/distance_template.h"
              *( (*(network + target)).buff_input_E + insert_position ) += ge;

              #include "templates/release_probability_template_stop.h"

            }

            tmp = ( (*(network+ii)).K + (int) floor( (double) K / (double) 3 ) );
            for ( jj = 0 ; jj < ( K % 3 ) ; jj++ )
            {

              #include "templates/release_probability_template_start.h"

              target = (uint32_t) ( ( tmp->offset & mask[jj] ) >> bitshift[jj] );
              #include "templates/distance_template.h"
              *( (*(network + target)).buff_input_E + insert_position ) += ge;

              #include "templates/release_probability_template_stop.h"

            }            

          }
          else
          {

            for ( jj = 0 ; jj < floor( (double) K / (double) 3 ) ; jj ++ )
            {

              tmp = ( (*(network+ii)).K + jj );

              #include "templates/release_probability_template_start.h"

              target = (uint32_t) ( ( tmp->offset & mask[0] ) >> bitshift[0] );
              #include "templates/distance_template.h"
              *( (*(network + target)).buff_input_I + insert_position ) += gi;

              #include "templates/release_probability_template_stop.h"

              #include "templates/release_probability_template_start.h"

              target = (uint32_t) ( ( tmp->offset & mask[1] ) >> bitshift[1] );
              #include "templates/distance_template.h"
              *( (*(network + target)).buff_input_I + insert_position ) += gi;

              #include "templates/release_probability_template_stop.h"

              #include "templates/release_probability_template_start.h"

              target = (uint32_t) ( ( tmp->offset & mask[2] ) >> bitshift[2] );
              #include "templates/distance_template.h"
              *( (*(network + target)).buff_input_I + insert_position ) += gi;

              #include "templates/release_probability_template_stop.h"

            }

            tmp = ( (*(network+ii)).K + (int) floor( (double) K / (double) 3 ) );
            for ( jj = 0 ; jj < ( K % 3 ) ; jj++ )
            {

              #include "templates/release_probability_template_start.h"

              target = (uint32_t) ( ( tmp->offset & mask[jj] ) >> bitshift[jj] );
              #include "templates/distance_template.h"
              *( (*(network + target)).buff_input_I + insert_position ) += gi;

              #include "templates/release_probability_template_stop.h"

            }

          }

        }

      }
      else
      {

        /* special for recording ^& */
        *( record_1_timestep + ii ) = (*(network+ii)).v;
        *( record_1_timestep_ge + ii ) = (*(network+ii)).ge;
        *( record_1_timestep_gi + ii ) = (*(network+ii)).gi;

      }

    }

    /* save spikes outside of record times too ^& */
    spike_array = start_spike_array;
    fwrite( spike_array, sizeof(int) + sizeof(double), spikes_in_this_round, o->spikefile );

    /* write values to file */
    /* special for recording ^& */
    if (  (tt >= record_start_bin) && (tt < record_stop_bin) &&
          (( (record_downsample_factor) == 0) || ((tt % record_downsample_factor) == 0 ))  )
    {

      /* write binned state variables to file */
  		if ( p->network.spatial_dimensions == 1 )
      {

    		for ( ii = 0 ; ii < vsteps ; ii++ )
    		{

     			sum_v = sum_ge = sum_gi = 0;
      		for ( jj = 0 ; jj < bin_size ; jj++ )
      		{
        		sum_v += *(record_1_timestep+(ii*bin_size)+jj);
        		sum_ge += *(record_1_timestep_ge+(ii*bin_size)+jj);
        		sum_gi += *(record_1_timestep_gi+(ii*bin_size)+jj);
      		}

        	*(record_1_timestep_binned+ii) = sum_v / bin_size;
    			*(record_1_timestep_binned_ge+ii) = sum_ge / bin_size;
    			*(record_1_timestep_binned_gi+ii) = sum_gi / bin_size;

    			totalsum_v += sum_v;
    			totalcount_v += bin_size;

        }

    	}
    	else if ( p->network.spatial_dimensions == 2 )
    	{

    		for ( ii = 0 ; ii < vsteps ; ii++ )
    		{

     			sum_v = sum_ge = sum_gi = 0;
      		origin = floor( ii / nvsteps ) * ( bin_size * NrowE ) + ( ii % nvsteps ) * ( bin_size );            
      		
      		for ( jj = 0 ; jj < bin_size ; jj++ )
          {
     	  		for ( kk = 0 ; kk < bin_size ; kk++ )
     	  		{
	        		sum_v += *(record_1_timestep+origin+(jj*NrowE)+kk);
	        		sum_ge += *(record_1_timestep_ge+origin+(jj*NrowE)+kk);
	        		sum_gi += *(record_1_timestep_gi+origin+(jj*NrowE)+kk);
    	  		}
      		}

      		*(record_1_timestep_binned+ii) = sum_v / (bin_size*bin_size);
          *(record_1_timestep_binned_ge+ii) = sum_ge / (bin_size*bin_size);
          *(record_1_timestep_binned_gi+ii) = sum_gi / (bin_size*bin_size);

          totalsum_v += sum_v;
          totalcount_v += (bin_size*bin_size);

      	}

    	}

      fwrite( record_1_timestep_binned, sizeof(double), vsteps, o->vfile );
      fwrite( record_1_timestep_binned_ge, sizeof(double), vsteps, o->gefile );
      fwrite( record_1_timestep_binned_gi, sizeof(double), vsteps, o->gifile );

    } 

    /* keep track of a few unbinned neurons to verify that everything looks normal */
    for ( jj = 0 ; jj < N ; jj+=record_individuals_jump )
    {
      fwrite( ( record_1_timestep+jj ), sizeof(double), 1, o->individual_traces );
    }
    /* end ^& */

  }

  /* write last values to the output struct */
  o->endtime = sim.T;
  o->vm_mean_all = totalsum_v / totalcount_v;

  /* return value */
  return NETSIM_NOERROR;

}

/* main function */
int main( int argc, char ** argv )
{

  /* print startup dialog */
  assert( startup_dialog() >= NETSIM_NOERROR );

  /* init */
  int ii;
  char * file_path = (char *)malloc( sizeof(char) * FILE_BUFFER_SIZE );
  clock_t begin, end;
  double buildtime, runtime, totaltime;
  struct parameters p = {{0}};

  /* allocate space for output path before parsing args */
  p.simulation.output_path = (char *)malloc( sizeof(char) * FILE_BUFFER_SIZE );
    
  /* parse input file path */
  p.simulation.parameter_file_path = (char *)malloc( sizeof(char) * FILE_BUFFER_SIZE );
  parse_inputs( argc, argv, &p );

  /* parse input file */
  p.simulation.optional_connection_path = (char *)malloc( sizeof(char) * FILE_BUFFER_SIZE );
  assert( parse_input_file( &p ) >= NETSIM_NOERROR );

  /* set defaults (not network parameters!) */
  if ( p.simulation.report_minutes == 0 ) {  p.simulation.report_minutes = 5;  }
  if ( p.neuron.vr != 0 && p.neuron.vreset == 0 ) {  p.neuron.vreset = p.neuron.vr;  }
  if ( p.network.poisson_rate>0 && p.network.poisson_stop_time == 0 ) {  p.network.poisson_stop_time = p.simulation.T;  }
  if ( p.network.poisson_rate == 0 ) {  p.network.poisson_stop_time = 0;  }
  if ( p.network.poisEIratio < .001 ) {  p.network.poisEIratio = .8;  }
  if ( p.network.spatial_dimensions == 0 ) {  p.network.spatial_dimensions = 1;  } /* 1D by default */
  if ( p.simulation.bin_size == 0 ) {  p.simulation.bin_size = 100;  }
  if ( p.simulation.cutoff_frequency == 0 ) {  p.simulation.cutoff_frequency = 160;  }
  if ( p.simulation.job_id == 0 ) {  p.simulation.job_id = INT_MIN;  }
  assert( (p.synapse.min_uniform_delay == 0) && (p.synapse.max_uniform_delay == 0) ); /* wrong simulator ^& */

  /* init network */
  struct neuron * network;
  network = malloc( sizeof(struct neuron) * p.network.N );

  /* init RNG */
  unsigned int seed = p.simulation.output_file_code;
  struct rng_state * rng;
  rng = malloc( sizeof(struct rng_state) );
  rng_init( rng, seed );

  /* init network */
  struct output o;
  o.spikecount = 0;
  assert( initialize_simulation( &p, network, &o, rng ) >= NETSIM_NOERROR );

  /* make connections */
  printf( "     Connecting network...\n" ); fflush( stdout ); begin = clock();

  if ( p.simulation.connector == NETSIM_RANDOM_CONNECT )  {  random_connect( &(p.network), network, rng );  }
  else if ( p.simulation.connector == NETSIM_GAUSSIAN_CONNECT )  {  gaussian_connect( &(p.network), &(p.simulation), &(p.synapse), network, rng );  }
  else if ( p.simulation.connector == NETSIM_TEST_CONNECT ) {  test_connect( &p, network );  }
  else if ( p.simulation.connector == NETSIM_GAUSSIAN_2D_CONNECT ) 
  {  
  	initialize_2D_simulation( &p, network ); 
    gaussian_connect_2D( &(p.network), &(p.simulation), &(p.synapse), network, rng );
    p.network.spatial_dimensions = 2;
  }  
  else {  printf("Given connector not specified\n"); return NETSIM_BAD_PARAMETER;  }

  /* random rewiring if specified */
  if ( p.network.rewiring_probability > 0 ) 
  {   
    printf( "     Rewiring network with probability %2.2f...\n", p.network.rewiring_probability ); fflush( stdout ); 
    random_rewiring( &p, network, rng ); 
  }

  /* save connectivity */
  if ( p.simulation.save_connectivity == 1 ) 
  	{  save_connectivity( &p, network ); save_incoming_connections( &p, network );  }

  end = clock(); buildtime = (double) (end - begin) / CLOCKS_PER_SEC;

  /* search for the right conductances to get desired rate */
  begin = clock();

  /* initialize input buffers */
  allocate_input_buffers( &p, network );

  /* is it a self-sustained simulation */
  if ( p.simulation.vm_mean != 0 || p.simulation.initiator == NETSIM_INIT_SUSTAINED )
  {  assert( initialize_sustained_simulation( &p, network, &o, rng ) <= NETSIM_NOERROR );  }

  /* start clock */
  printf( "Running network with ge %.10f nS and gi %.10f nS\n", p.synapse.ge*1e9, p.synapse.gi*1e9 );
  fflush( stdout );
  begin = clock();

  /* run simulation */
  assert( run_simulation( &p, network, &o, rng ) <= NETSIM_NOERROR );

  /* calculate simulation time */
  end = clock(); runtime = (double)( end - begin ) / CLOCKS_PER_SEC;
  totaltime = buildtime + runtime;

  /* print summary values to screen */
  printf( "Endtime: %g\nRate (Hz): %g\nBuildtime (s): %g, Runtime (s): %g, Totaltime (s): %g\n",
    o.endtime, ((double)o.spikecount)/p.network.N/o.endtime, buildtime, runtime, totaltime );
  printf( "Total number of spikes: %d\n", o.spikecount );
  printf( "Mean membrane potential: %f (mV)\n", (o.vm_mean_all)*1000 );

  /* clean up */
  fclose( o.spikefile );
  fclose( o.vfile );
  fclose( o.gefile );
  fclose( o.gifile );
  fclose( o.individual_traces );
  fclose( o.individual_traces_ge );

  for ( ii = 0 ; ii < p.network.N ; ii++ ) {  free( (*(network+ii)).K );  }
  free( network ); free( rng );
  free( file_path );

  free( p.simulation.optional_connection_path );
  free( p.simulation.output_path );
  free( p.simulation.parameter_file_path );

  return NETSIM_NOERROR;

}
