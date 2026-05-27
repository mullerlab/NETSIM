#include <assert.h>
#include <getopt.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h> // Include for mkdir

#include "isaac64.h"
#include "rng.h"
#include "rng_dist.h"

#include "connectors.h"
#include "initializers.h"
#include "input.h"
#include "netsim.h"
#include "output.h"

/* define constants */
#define NETSIM_DIALOG_NOCOLOR__ "\x1b[0m"
#define NETSIM_DIALOG_INVCOLOR__ "\x1b[7m"

#define NETSIM_DIALOG_BG_COLOR__RED__ "\x1b[41;1m"

#define NETSIM_DIALOG_FG_COLOR__RED__ "\x1b[31;1m"
#define NETSIM_DIALOG_FG_COLOR__GREEN__ "\x1b[32;1m"
#define NETSIM_DIALOG_FG_COLOR__YELLOW__ "\x1b[33;1m"
#define NETSIM_DIALOG_FG_COLOR__BLUE__ "\x1b[34;1m"

/* NETSIM startup dialog */
int startup_dialog(void) {

    printf("\n\n");

    printf("===============================\n");
    printf("============ %sNETSIM%s ===========\n", NETSIM_DIALOG_FG_COLOR__RED__,
           NETSIM_DIALOG_NOCOLOR__);
    printf("===============================\n");

    printf("\n\n");

    return NETSIM_NOERROR;
}

#include <sys/stat.h> // Include for mkdir

int parse_inputs(int argc, char **argv, struct parameters *p) {
    /* init options */
    int option = 0;

    /* parse inputs */
    while ((option = getopt(argc, argv, "f:j:o:")) != -1) {
        switch (option) {
        case 'f':
            strcpy(p->simulation.parameter_file_path, optarg);
            break;
        case 'j':
            p->simulation.job_id = (unsigned int)atoi(optarg);
            break;
        case 'o':
            strcpy(p->simulation.output_path, optarg);
            p->simulation.output_path_override = 1;
            printf("Command line output_path: %s overriding parameter file\n", optarg);

            // Check if the directory exists, and create it if it doesn't
            struct stat st = {0};
            if (stat(optarg, &st) == -1) {
                int result = mkdir(optarg, 0700);
                assert(result == 0 && "Error: Could not create output directory");
                printf("Created output directory: %s\n", optarg);
            }
            break;
        default:
            printf("USAGE: \"netsim -j job_id -f parameter_file -o output_directory\"\n");
            return NETSIM_BAD_PARAMETER;
        }
    }

    /* return */
    return NETSIM_NOERROR;
}

int param_error_checking(struct parameters p) {
    bool current_injection = false;
    bool external_input = false;
    bool variable_external_input = false;
    bool release_probability = false;
    bool long_simulation = true;

#if defined(CURRENT_INJECTION)
    current_injection = true;
#endif

#if defined(EXTERNAL_INPUT)
    external_input = true;
#endif

#if defined(VARIABLE_EXTERNAL_INPUT)
    variable_external_input = true;
#endif

#if defined(RELEASE_PROBABILITY)
    release_probability = true;
#endif

#if defined(LONG_SIMULATION)
    long_simulation = true;
#endif
    
    // external input
    if (p.network.poisson_start_time < p.network.poisson_stop_time && p.network.poisson_rate != 0) {
        assert(external_input);
    }

    // variable external input
    if (p.network.poisson_start_time < p.network.poisson_stop_time &&
        strcmp(p.simulation.poisson_rate_path, "") != 0) {
        assert(variable_external_input);
    }

    // need number of current sources and the path
    if (p.simulation.current_injection_count > 0) {
        assert(p.network.N % p.simulation.current_injection_count == 0);
        assert(current_injection);
    } else {
        assert(!current_injection);
    }

    if(p.simulation.partial_recording_time != 0) {
	    assert(long_simulation);
    }
    
    /* set defaults (not network parameters!) */
    if (p.neuron.vr != 0 && p.neuron.vreset == 0) {
        p.neuron.vreset = p.neuron.vr;
    }
    if (p.network.spatial_dimensions == 0) {
        p.network.spatial_dimensions = 1;
    } /* 1D by default */
    if (p.simulation.job_id == 0) {
        p.simulation.job_id = INT_MIN;
    }
    assert((p.synapse.min_uniform_delay == 0) &&
           (p.synapse.max_uniform_delay == 0)); /* wrong simulator ^& */
    
    assert(p.simulation.start_record_time >= 0.0 &&
 	   p.simulation.start_record_time < p.simulation.T &&
 	   p.simulation.stop_record_time >= 0.0 &&
	   p.simulation.stop_record_time <= p.simulation.T &&
 	   p.simulation.start_record_time <= p.simulation.stop_record_time);
    
    return 0;
}

/* main function */
int main(int argc, char **argv) {
    /* print startup dialog */
    assert(startup_dialog() >= NETSIM_NOERROR);

    /* init */
    struct parameters p = {{0}};

    /* allocate space for output path before parsing args */
    p.simulation.output_path = calloc(FILE_BUFFER_SIZE, sizeof(char));
    p.simulation.current_injection_path = calloc(FILE_BUFFER_SIZE, sizeof(char));
    p.simulation.sparse_current_injection_path = calloc(FILE_BUFFER_SIZE, sizeof(char));
    p.simulation.poisson_rate_path = calloc(FILE_BUFFER_SIZE, sizeof(char));
    /* make sure they are empty strings */
    *p.simulation.output_path = '\0';
    *p.simulation.current_injection_path = '\0';
    *p.simulation.sparse_current_injection_path = '\0';
    *p.simulation.poisson_rate_path = '\0';

    /* parse input file path */
    char *file_path = calloc(FILE_BUFFER_SIZE, sizeof(char));
    p.simulation.parameter_file_path = calloc(FILE_BUFFER_SIZE, sizeof(char));
    *p.simulation.parameter_file_path = '\0';
    assert(parse_inputs(argc, argv, &p) == NETSIM_NOERROR);

    /* parse input file */
    p.simulation.optional_connection_path = calloc(FILE_BUFFER_SIZE, sizeof(char));
    *p.simulation.optional_connection_path = '\0';
    assert(parse_input_file(&p) == NETSIM_NOERROR);

    param_error_checking(p);

    /* init network */
    struct neuron *network = NULL;
    network = calloc(1, sizeof(struct neuron));

    /* init RNGs */
    unsigned int seed = p.simulation.output_file_code;
    struct rng_state *rng = NULL;
    rng = calloc(1, sizeof(struct rng_state));
    rng_init(rng, seed);

    #ifdef USE_SEPARATE_CONNECTOR_RNG
        unsigned int conn_seed = p.simulation.connection_seed;
        struct rng_state *conn_rng = NULL;
        conn_rng = malloc(sizeof(struct rng_state));
        rng_init(conn_rng, conn_seed);
    #endif

    clock_t begin = clock();

    /* init network */
    struct output o;
    o.spike_count = 0;

    /* make connections */
    assert(initialize_network(&p, network, rng) == NETSIM_NOERROR);

    #ifdef USE_SEPARATE_CONNECTOR_RNG
        assert(connect_network(&p, network, conn_rng) == NETSIM_NOERROR);
    #else
        assert(connect_network(&p, network, rng) == NETSIM_NOERROR);
    #endif

    /* init external inputs*/
    assert(init_simulated_inputs(&p) == NETSIM_NOERROR);

    /* save connectivity */
    if (p.simulation.save_connectivity == 1) {
        save_connectivity(&p, network);
    }

    /* initialize input buffers */
    allocate_input_buffers(&p, network);

    clock_t end = clock();
    double buildtime = (double)(end - begin) / CLOCKS_PER_SEC;

    begin = clock();

    /* is it a self-sustained simulation */
    if (p.simulation.vm_mean != 0 || p.simulation.initiator == NETSIM_INIT_SUSTAINED) {
        assert(initialize_sustained_simulation(&p, network, rng) == NETSIM_NOERROR);
    }

    /* start clock */
    printf("Running network with ge %.10f nS and gi %.10f nS\n", p.synapse.ge * 1e9,
           p.synapse.gi * 1e9);
    fflush(stdout);
    begin = clock();

    /* run simulation */
    assert(run_simulation(&p, network, &o, rng) == NETSIM_NOERROR);

    /* calculate simulation time */
    end = clock();
    double runtime = (double)(end - begin) / CLOCKS_PER_SEC;
    double totaltime = buildtime + runtime;

    /* print summary values to screen */
    printf("Endtime: %g\nRate (Hz): %g\nBuildtime (s): %g, Runtime (s): %g, Totaltime (s): %g\n",
           p.simulation.T,
	   ((double)o.spike_count)/p.network.N/p.simulation.T,
	   buildtime,
	   runtime,
           totaltime);
    printf("Total number of spikes: %zu\n", o.spike_count);

    /* save out parameters as they exist */
    save_out_param_file(&p, o.spike_count);

    // free circular buffers
    free(network->buff_input_E);
    free(network->buff_input_I);

    // free neuron structs
    free(network->excitatory);
    free(network->v);
    free(network->ge);
    free(network->gi);
    free(network->K);
    free(network->conduction_speed);
    free(network->x_position);
    free(network->y_position);
    free(network->z_position);
    free(network->lastSpike);

    //
    free(network);
    free(rng);
#ifdef USE_SEPARATE_CONNECTOR_RNG
    free(conn_rng);
#endif
    free(file_path);

    //
    free(p.simulation.optional_connection_path);
    free(p.simulation.output_path);
    free(p.simulation.parameter_file_path);
    free(p.simulation.current_injection_path);
    free(p.simulation.sparse_current_injection_path);
    free(p.simulation.poisson_rate_path);

    // if used, free current injection structs
    if (p.simulation.sparse_current_entry_count > 0) {
        for (int ii = 0; ii < p.simulation.sparse_current_entry_count; ii++) {
            free(p.simulation.sparse_injection[ii].indices);
        }
        free(p.simulation.sparse_injection);
    }

    if (p.simulation.current_injection_count > 0) {
        free(p.simulation.current_injection);
    }

    /* return */
    return NETSIM_NOERROR;
}
