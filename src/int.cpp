/*	int.c
	routines for numerical integration
	started 1 June 90, es using old FORTRAN stuff and newer C stuff
*/

// This file contains all the integration routines needed.

#include	<stdio.h>
#include	<math.h>
#include        <stdlib.h>
#include	"int.h"

#ifdef TEST


Int::main() {
	testgala();
	}

int	Int::matherr() {		/* for debugging */
	printf("matherr\n");	
	}
	
Int::testgauss() {
	double	xlo = 0, xhi = 0.99;
	printf( "quadrature gives exact %10f:\n", asin(xhi) - asin(xlo) );
	printf( "%2dpoints: %10f\n", 4,gauss( 4,testfunc,xlo,xhi) );
	printf( "%2dpoints: %10f\n", 8,gauss( 8,testfunc,xlo,xhi) );
	printf( "%2dpoints: %10f\n",10,gauss(10,testfunc,xlo,xhi) );
	printf( "%2dpoints: %10f\n",12,gauss(12,testfunc,xlo,xhi) );
	printf( "%2dpoints: %10f\n",16,gauss(16,testfunc,xlo,xhi) );
	printf( "%2dpoints: %10f\n",20,gauss(20,testfunc,xlo,xhi) );
	printf( "%2dpoints: %10f\n",48,gauss(48,testfunc,xlo,xhi) );
	}

Int::testgala() {
	printf( "testing Gauss-Laguerre integration with sin(2x)*exp(-2x)\n" );
	printf( "quadrature gives exact %10f:\n", 0.25 );
	printf( "%2dpoints: %10f\n", 4, gala( 4,testgalafunc,0.0,0.5) ); 
	printf( "%2dpoints: %10f\n", 8, gala( 8,testgalafunc,0.0,0.5) );
	printf( "%2dpoints: %10f\n",12, gala(12,testgalafunc,0.0,0.5) );
	}

Int::testgahe() {
	printf( "testing Gauss-Hermite integration with exp(-x*x)\n" );
	printf( "exact result: %10f\n", sqrt( M_PI ) );
	printf( "%2d points: %10f\n", 4, gahe( 4,testfunc,0.0,1.0) );
	printf( "%2d points: %10f\n", 8, gahe( 8,testfunc,0.0,1.0) );
	printf( "%2d points: %10f\n",16, gahe(16,testfunc,0.0,1.0) );
	}

Int::testgauche() {
	double	xpole = 100, xbase = 0;
	printf("%.15f\n", asin(1.0) * 2 );
	printf( "testing Gauss-Chebyshev integration with 1/sqrt(1-x*x)\n" );
	printf( "exact result:%10f\n", asin(xpole)-asin(xbase) );
	printf( "  %2d points: %10f\n", 4, gauche( 4,testfunc,xbase,xpole) );
	printf( "  %2d points: %10f\n", 8, gauche( 8,testfunc,xbase,xpole) );
	printf( "  %2d points: %10f\n",16, gauche(16,testfunc,xbase,xpole) );
	printf( "  %2d points: %10f\n",96, gauche(96,testfunc,xbase,xpole) );
	}

double	Int::testgalafunc( x )
	double	x;
	{
	return ( sin(2*x)*exp(-2*x) );
	}

double Int::testfunc( x )
	double x;
	{
	/* x += 0.5; return exp( -x*x );	*/
	return 1/sqrt(100*100-x*x);
	}
#endif



/******************************************************
*
*	gauss
*
*
* Gauss-Legendre Quadrature w/ switchable no of points 
* 4 Jun 90, es
********************************************************/

double Int::gauss(int n, double (*f)(double, void *), double xlo, double xhi, void *optvec )
	{
	double	xoffs, xdiff; 
	int	ix;
	double	s;		/* summing up */
	double	*p, *w;		/* pointing to active list */

	switch ( n ) {
		case 4:		p= gaulep4; w= gaulew4; break;
		case 8:		p= gaulep8; w= gaulew8; break;
		case 10:	p=gaulep10; w=gaulew10; break;
		case 12:	p=gaulep12; w=gaulew12; break;
		case 16:	p=gaulep16; w=gaulew16; break;
		case 20:	p=gaulep20; w=gaulew20; break;
		case 48:	p=gaulep48; w=gaulew48; break;
		default:	printf("\ngauss():%d points not in list\n",n);
				exit(0);
		}
	xoffs = 0.5 * ( xlo + xhi );
	xdiff = 0.5 * ( xhi - xlo );
	s = 0;
	for( ix=0; ix<n/2; ix++ ) 	/* n is even */
		s += w[ix] * ( f(xoffs+xdiff*p[ix],optvec)
			     + f(xoffs-xdiff*p[ix],optvec) );
	return( s * xdiff );
	}


/******************************************************
*
*	gausspts
*
*
* returns points for Gauss-Legendre Quadrature
* from gauss(), 25 May 91, es
********************************************************/

void	Int::gausspts( n, xlo, xhi, xvec, wvec )
	int	n;		/* number of points must be even */
	double	xlo, xhi;	/* limits */
	double	*xvec, *wvec;	/* abszissas and weights	*/
	{
	double	xoffs, xdiff; 
	int	ix;
	double	s;		/* summing up */
	double	*p, *w;		/* pointing to active list */

	switch ( n ) {
		case 4:		p= gaulep4; w= gaulew4; break;
		case 8:		p= gaulep8; w= gaulew8; break;
		case 10:	p=gaulep10; w=gaulew10; break;
		case 12:	p=gaulep12; w=gaulew12; break;
		case 16:	p=gaulep16; w=gaulew16; break;
		case 20:	p=gaulep20; w=gaulew20; break;
		case 48:	p=gaulep48; w=gaulew48; break;
		default:	printf("\ngauss():%d points not in list\n",n);
				exit(0);
		}
	xoffs = 0.5 * ( xlo + xhi );
	xdiff = 0.5 * ( xhi - xlo );
	for( ix=0; ix<n/2; ix++ ) {	/* n is even */
		xvec[ix]	= xoffs-xdiff*p[ix];
		xvec[n-1-ix]	= xoffs+xdiff*p[ix];
		wvec[ix]	= wvec[n-1-ix] = xdiff * w[ix];
		}
	}

/******************************************************
*
*	gaussn
*
*
* Gauss-Legendre Quadrature w/switchable no of points 
* 	subdivisions of region, parameter handover
* 11 Sep 90, es
********************************************************/

double	Int::gaussn( n, ndiv, f, xlo, xhi, optvec )
	int	n;		/* number of points must be even	*/
	int	ndiv;		/* number of subdivisions		*/
	double	(*f)();		/* function of one double parameter	*/
	double	xlo, xhi;	/* limits				*/
	void	*optvec;	/* parameters, hand over to function	*/
	{
	double	xoffs, xdiff; 
	int	ix, idiv;
	double	s;		/* summing up */
	double	*p, *w;		/* pointing to active list */

	switch ( n ) {
		case 4:		p= gaulep4; w=gaulew4; break;
		case 8:		p= gaulep8; w=gaulew8; break;
		case 10:	p=gaulep10; w=gaulew10; break;
		case 12:	p=gaulep12; w=gaulew12; break;
		case 16:	p=gaulep16; w=gaulew16; break;
		case 20:	p=gaulep20; w=gaulew20; break;
		case 48:	p=gaulep48; w=gaulew48; break;
		default:	printf("\ngauss():%d points not in list\n",n);
				exit(0);
		}
	s = 0;
	xdiff = 0.5 * ( xhi - xlo ) / ndiv;
	for( idiv=1; idiv<=ndiv; idiv++ ) {
		xoffs = xlo + (2*idiv-1) * xdiff;
		for( ix=0; ix<n/2; ix++ ) 	/* n is even */
			s += w[ix] * ( f(xoffs+xdiff*p[ix],optvec) 
				     + f(xoffs-xdiff*p[ix],optvec) );
		}
	return( s * xdiff );
	}

/***************************************************
*
*	gala()
*
* Gauss-Laguerre quadrature extends to +oo
* adapted from old FORTRAN, 24 JUL 90, es, 17 aug 90, 27 Mar 91
****************************************************/
double	Int::gala( n, f, xlo, invslope, optvec )
	int	n;	/* number of points		*/
	double	(*f)();	/* function to integrate	*/
	double	xlo;	/* lower limit of integration	*/
	double	invslope;/* approximate inverse slope of decaying 
			   function is roughly proportional to the
			   integration region, needed for placing 
			   the points. i.e. =1 for exp(-x), 
			   invslope=0.5 for exp(-2), etc	*/
	void	*optvec; /* optional vector, passesd to function f() */
	{
	double	*x, *w;
	int	i;
	double	sum	=0;

	if( n == 4 )		{ x = gala4x; w = gala4w; }
	else if( n == 8 )	{ x = gala8x; w = gala8w; }
	else if( n == 12 )	{ x = gala12x;w = gala12w;}
	else if( n == 15 )	{ x = gala15x;w = gala15w;}
	else {
		printf("\ngala():n=%d not in list\n", n );
		return 0;
		}
	for( i=0; i<n; i++ ) 
		sum += w[i] * f( invslope*x[i] + xlo, optvec );
	return invslope * sum;	/* make up for transformation */
	}

/******************************************************
*
*	gahe
*
*
* Gauss-Hermite Quadrature w/ switchable no of points 
* 19 Sep 90, es
********************************************************/
double	Int::gahe( n, f, center, width, optvec )
	int	n;		/* number of points must be even */
	double	(*f)();		/* function of one double parameter */
	double	center;		/* approx center of integr. region */
	double	width;		/* approx width for integration region */
	void	*optvec;	/* optional vector, passed to function	*/
	{
	int	ix;
	double	dx;
	double	s;		/* summing up */
	double	*p, *w;		/* pointing to active list */

	switch ( n ) {
		case 4:		p= gahep4;  w= gahew4; break;
		case 8:		p= gahep8;  w= gahew8; break;
		case 16:	p= gahep16; w=gahew16; break;
		case 20:	p= gahep20; w=gahew20; break;
		default:	fprintf(stderr,
					"\ngahe():%d points not in list\n",n);
				exit(0);
		}
	s = 0;
	for( ix=0; ix<n/2; ix++ ) {	/* n is even */
		dx = width * p[ix];
		s += w[ix] * ( f(center+dx,optvec) + f(center-dx,optvec) );
		}
	return s * width;
	}


/************************************************************************
*
*	gauche
*
*
* Gauss-Chebyshev Quadrature w/no.points and parameters 
* weight function is (1-x*x)^{-1/2} with x in [0,1] intervall only
* 27 Sep 90, es (see Abramowitz 25.4.38)
*************************************************************************/
double	Int::gauche( n, f, base, pole, optvec )
	int	n;			/* number of points (is free)	*/
	double	(*f)();			/* function to be integrated	*/
	double	base;			/* one endpoint of region	*/
	double	pole;			/* other endpoint is a pole	*/
	void	*optvec;		/* just passed on to function	*/
					/* is optional argument		*/
	{	
	static	int	list[16] = {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
					/* which list has been computed	*/
	static	double	*plist[16];	/* the points where to evaluate	*/
	static	double	*wlist[16];	/* the weights there		*/
	double	*p, *w;			/* points/weights for work	*/
	double	dist	= pole - base;	/* need not be positive		*/
	double	sum = 0, sumsum = 0;	/* summing in bundels		*/
	int	i;

	i = 0;
	while( i<16 && list[i] != 0  && list[i] != n )
		i++;
	if( i==16 ) {
		fprintf( stderr, "\ngauche():list full (max 16)\n");
		exit( 1 );
		}
	if( list[i] == 0 ) {	/* not in list yet, construct */
		p = ((double *) malloc( (unsigned)n*sizeof(double) ));
		w = ((double *) malloc( (unsigned)n*sizeof(double) ));
		if( p && w )		/* alloc worked */
			{ p--; w--; }	/* make 1-offset vectors */
		else {
			fprintf( stderr, "\ngauche():malloc failed\n" );
			exit( 1 );
			}
		list[i]  = n;
		plist[i] = p;
		wlist[i] = w;
		for( i=1; i<=n; i++ ) {		/* take one half only	*/
			p[i] = cos( (2*i-1)*M_PI/(4*n) );
			w[i] = sin( (2*i-1)*M_PI/(4*n) ); /* M_PI/(2*n) */
			}
		}
	else {		/* already in list */
		p = plist[i];
		w = wlist[i];
		}
	for( i=1; i<=n; i++ ) {
		sum += w[i] * f( base + dist * p[i], optvec );
		if( i&0xfff0 == 0 ) {		/* sum in bundels	*/
			sumsum += sum;
			sum = 0;
			}
		}
	sumsum += sum;
	return	dist * sumsum * M_PI / (2*n); /* rest of weight factor */
	}
		

/******************************************************
*
*	gaussp, galap, gahep
*
* just kept for compatibility, 27 Mar 91, es
********************************************************/

double	Int::gaussp( n, f, xlo, xhi, para )
	int	n;
	double	(*f)();	
	double	xlo, xhi;
	double	para[];	
	{
	return	gauss( n, f, xlo, xhi, para );
	}

double	Int::galap( n, f, xlo, invslope, para )
	int	n;
	double	(*f)();
	double	xlo;
	double	invslope;
	double	para[];
	{
	return gala( n, f, xlo, invslope, para );
	}

double	Int::gahep( n, f, center, width, para )
	int	n;		/* number of points must be even */
	double	(*f)();		/* function of one double parameter */
	double	center;		/* approx center of integr. region */
	double	width;		/* approx width for integration region */
	double	*para;		/* parameter block */
	{
	return gahe( n, f, center, width, para );
	}


/******************************************
*
*	gaussnbyn
*
*
* twodimensional Gauss-Legendre Quadrature
* 15 Mar 90, es
*******************************************/

double	Int::gaussnbyn( n, f, xlo, xhi, ylo, yhi )
	int	n;	/* number of points per direction */
	double	(*f)();	/* function of two parameters	*/
	double	xlo, xhi, ylo, yhi;	/* limits */
	{
	double	xoffs, xdiff, x; 
	double	yoffs, ydiff, y;
	int	ix, iy;
	double	wx, sy, s;
	double	*p, *w;		/* pointing to active list */
	static	double	p10list[] = {	0.1488743389,	0.4333953941, 
			0.6794095682,	0.8650633666,	0.97390652	};
	static	double	w10list[] = {	0.2955242247,	0.2692667193,
			0.2190863625,	0.1494513491,	0.06667134	};

	if( n==10) {
		p = p10list;
		w = w10list;
		}
	else  
		printf( "\ngaussnbyn(): this number not in list\n" );

	xoffs = 0.5 * ( xlo + xhi );
	yoffs = 0.5 * ( ylo + yhi );
	xdiff = 0.5 * ( xhi - xlo );
	ydiff = 0.5 * ( yhi - ylo );
	s = 0;
	for( ix=0; ix<n; ix++ ) {
		if( ix<n/2 ) {
			x = xoffs + xdiff * p[ix];
			wx = w[ix];
			}
		else {
			x = xoffs - xdiff * p[n-1-ix];
			wx = w[n-1-ix];
			}
		sy = 0;
		for( iy=0; iy<n/2; iy++ ) 
			sy += w[iy] * 
			( f(x,yoffs+ydiff*p[iy]) + f(x,yoffs-ydiff*p[iy]) );
		s += wx * sy;
		}
	return( s * xdiff * ydiff );
	}
