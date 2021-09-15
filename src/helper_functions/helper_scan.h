/*

helper functions - parameter scans

*/

double linspace_val( double start, double stop, int npts, int val_index )
{

	/* init */
	double val, interval;

	/* create output */
	interval = ( stop - start ) / ( npts - 1 );
	val = start + (double)(val_index)*interval;

	/* return */
	return val;

}

int sub2ind( unsigned int x1, unsigned int x2, unsigned int x3, unsigned int n1, unsigned int n2, unsigned int n3 )
{

	/* init */
	unsigned int index = 0;

	/* error checking */
	assert( (x1<=n1) && (x2<=n2) && (x3<=n3) );

	/* calculate linear index output */
	index = x1 + (x2-1)*n2 + (x3-1)*n1*n2;

	/* return */
	return index;

}

int ind2sub( unsigned int index, unsigned int dim, unsigned int n1, unsigned int n2, unsigned int n3 )
{

	/* init */
	unsigned int subscript_dim = 0;

	/* error checking */
	assert( index <= ((n1+1)*(n2+1)*(n3+1)) );	

	/* calculate subscript output */
	if ( dim == 1 ) {  subscript_dim = index % n1;  }
	else if ( dim == 2 ) {  subscript_dim = floor( index / n1 );  }
	else if ( dim == 3 ) {  subscript_dim = floor( index / (n1*n2) ); }

	/* return */
	return subscript_dim;

}
