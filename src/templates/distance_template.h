
/* special for recalculating delay @# */
distance_x = fabs( ( *(network+ii) ).x_position - ( *(network + target ) ).x_position );
distance_x = fmin( L - distance_x, distance_x );
distance_y = fabs( ( *(network+ii) ).y_position - ( *(network + target ) ).y_position );
distance_y = fmin( L - distance_y, distance_y );

distance = sqrt( distance_x*distance_x + distance_y*distance_y );

delay_in_bins = (int) round( ( synapse_delay + ( distance / ( *(network+ii) ).conduction_speed ) ) / dt );
insert_position = ( buffer_position + delay_in_bins ) % buffer_length;
/* end @# */
