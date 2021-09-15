/*

helper functions - gaussian connect

*/

//// DEFINITELY PROFILE
double calculate_shortest( int i, int j, double dxs, double dxr, double L ) 
{
  // Calculates shortest distance bewteen two points on a ring with length L 
  // if distance bewteen each point is equally spaced, dx
  
  double opt1 = fabs(i*dxs - j*dxr);
  double opt2 = (L - opt1);
  return (( opt1 < opt2 ) ?  opt1 : opt2 );

}


/* function: calculate shortest distance between two nodes on a circle. */
/* warning: distance given has satisfy  distance/dxj > -N */
int calculate_closest( double distance, double position, double dxj, int N )
{

  double jdis;
  int j;

  jdis = position + distance;
  j = round ( jdis/dxj );
  
  j = (j+N)%N;

  return j;
  
}