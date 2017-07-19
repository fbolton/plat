#include "include.h"
#include <sys/time.h>

/*********************/
/* Hexagonal network */
/*********************/

#define hexv(i,j) ((i)*(jmax)+(j))
#define hexc(i,j) ((i)*(jmax/2)+((j)/2))

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *    H E X N E T ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	hexnet(short nx, short ny, REAL f)
 *
 *	Arguments:	nx	= number of hex cells vertically
 *			ny	= 2 * number of hex cells horizontally
 *			f	= fraction of cell edge covered by Plateau
 *				border (i.e. lies in range [0.0, 1.0] )
 *
 *	Return value:	none
 *
 *	Action:		Sets up a periodic network of regular hexagonal cells
 *			(i.e. a honeycomb) including Plateau borders.  The
 *			network consists of 2 * ny cells horizontally, before
 *			the period wraps around, and `nx' cells vertically.
 *			The total size of the network will thus be `2*nx*ny'.
 *			Plateau borders have their size determined by
 *			`f' which is the fraction of a cell edge covered.
 *
 *			The construction of the froth is done in two stages.
 *			Initially the network is constructed as a *dry* froth
 *			with no Plateau borders, and vertices placed at the
 *			apices of perfect hexagons.  These vertices are
 *			organised approximately on a 2-d grid which is then
 *			indexed by the integers `i' and `j'.  This 2-d
 *			indexing scheme is mapped to the 1-d index of the
 *			vertex arrays (e.g. vx[], vy[]) via the macro
 *			`hexv(i,j)' defined above.  The slightly messy
 *			procedure of defining the topology is then taken care
 *			of within a 2-d loop (e.g. the topological arrays
 *			`vorvnbr[][]', `vorvper[][]' and `vorcadj[][]' are
 *			defined.  Note that these arrays are specially used to
 *			store the topology of the dry network and are *not* the
 *			arrays ultimately used to store the Plateau border
 *			topology).
 *
 *			The second stage is accomplished via calling the
 *			general routine `vorplateautopol()' which converts a
 *			dry froth to a froth with Plateau borders.  The
 *			variables after this point refer to the Plateau border
 *			network topology.
 *
 *****************************************************************************/
void hexnet(nx,ny,f)
short nx,ny;
REAL f;  /* f = fraction of edge covered by PB */
{
  short nmax, imax, jmax, ii, k, c1, i, j, ij, ij1, id, jd, iu, ju, iuu, juu;
  short csides();
  REAL rt3, a, b, a0, b0, yshift, d, cp1, ca1, ba1, cellarea(), sqf();
  boolean found[MCELL], arcbrk;
  void vorplateautopol(), plerror();
  nx= (nx>0) ? nx : 1;
  ny= (ny>0) ? ny : 1;
  if (f<0 || f>1.0) { plerror("hexnet: f out of range"); return; }
  f= (f==0) ? 0.001 : f;
  rt3=sqrt(3.0);
  if ( (a=rt3*nx) > (b=((REAL) ny)) ) {
    boxwid=1.0; boxhgt=b/a;
  }
  else {
    boxwid=a/b; boxhgt=1.0;
  }
  a=boxwid/(2.0*nx);
  b=boxhgt/(2.0*ny);
  yshift= -0.1*b;
  imax=2*nx; jmax=2*ny;
  for (i=0; i<imax; i+=2)
    for (j=0; j<jmax; j+=2) {
      /* See defn. of hexv() above to see the relation */
      /* between the 1- and the 2-dim indexing. */
      a0=a*(i-nx+1); b0=b*(j-ny+1);
      id=(i-1+imax)%imax; jd=(j-1+jmax)%jmax;
      iu=(i+1)%imax; ju=(j+1)%jmax;
      iuu=(i+2)%imax; juu=(j+2)%jmax;
      vx[hexv(i,ju)]=a0-a/3.0;
      vx[hexv(iu,ju)]=a0+a/3.0;
      vy[hexv(i,ju)]=vy[hexv(iu,ju)]=b0+b+yshift;
      vx[hexv(i,j)]=a0-2.0*a/3.0;
      vx[hexv(iu,j)]=a0+2.0*a/3.0;
      vy[hexv(i,j)]=vy[hexv(iu,j)]=b0+yshift;
      ij=hexv(i,j);
      vorvnbr[ij][0]=hexv(id,j);
      vorvper[ij][0]= (id<i) ? 0 : PERFN(-1,0);
      vorcadj[ij][0]=hexc(i,j);
      vorvnbr[ij][1]=hexv(i,jd);
      vorvper[ij][1]= (jd<j) ? 0 : PERFN(0,-1);
      vorcadj[ij][1]=hexc(id,j);
      vorvnbr[ij][2]=hexv(i,ju);
      vorvper[ij][2]=0;
      vorcadj[ij][2]=hexc(id,jd);
      ij=hexv(i,ju);
      vorvnbr[ij][0]=hexv(iu,ju);
      vorvper[ij][0]=0;
      vorcadj[ij][0]=hexc(id,j);
      vorvnbr[ij][1]=hexv(i,juu);
      vorvper[ij][1]= (juu>ju) ? 0 : PERFN(0,1);
      vorcadj[ij][1]=hexc(i,j);
      vorvnbr[ij][2]=hexv(i,j);
      vorvper[ij][2]=0;
      vorcadj[ij][2]=hexc(i,juu);
      ij=hexv(iu,ju);
      vorvnbr[ij][0]=hexv(i,ju);
      vorvper[ij][0]=0;
      vorcadj[ij][0]=hexc(iu,j);
      vorvnbr[ij][1]=hexv(iu,j);
      vorvper[ij][1]=0;
      vorcadj[ij][1]=hexc(i,juu);
      vorvnbr[ij][2]=hexv(iu,juu);
      vorvper[ij][2]= (juu>ju) ? 0 : PERFN(0,1);
      vorcadj[ij][2]=hexc(i,j);
      ij=hexv(iu,j);
      vorvnbr[ij][0]=hexv(iuu,j);
      vorvper[ij][0]= (iuu>iu) ? 0 : 1;
      vorcadj[ij][0]=hexc(i,j);
      vorvnbr[ij][1]=hexv(iu,ju);
      vorvper[ij][1]=0;
      vorcadj[ij][1]=hexc(iu,jd);
      vorvnbr[ij][2]=hexv(iu,jd);
      vorvper[ij][2]= (jd<j) ? 0 : PERFN(0,-1);
      vorcadj[ij][2]=hexc(iu,j);
    }
  nc=2*nx*ny; nv=2*nc;
  d=f*a/3.0;
  vorplateautopol(d);
  for (i=0; i<nv; i++) vlist[i]=i;
  for (i=0; i<nc; i++) clist[i]=i;
  for (i=0; i<nb; i++) blist[i]=i;
  onv=nv; onc=nc; onb=nb;
  bpav= -tan(PI/6.0)/d;
  for (i=0; i<nb; i++) bp[i]=bpav;
  ba1=(sqrt(3.0)-PI/2.0)*sqf(BRADIUS(0.0,bpav));
  for (i=0; i<nb; i++) barea[i]=ba1;
  for (i=0; i<nc; i++) {
    darea[i]=0.0; cp[i]=0.0; ncsides[i]=6;
  }
/*ca1=boxwid*boxhgt/nc;*/
/*for (i=0; i<nc; i++) carea[i]=ca1;*/
  for (i=0; i<onc; i++) found[i]=FALSE;
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    for (k=1; k<3; k++) {
      if (!found[c1=cadj[i][k]]) {
        found[c1]=TRUE;
        carea[c1]=cellarea(i,k,&arcbrk);
        if (arcbrk) { plerror("failure of routine hexnet"); return; }
      }
    }
  }
  henckyeps=0.0;
}

/**********************/
/* 'Voronoi' routines */
/**********************/

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     * V O R S P R A Y ( )  *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	vorspray(REAL *x, REAL *y, REAL wid, REAL hgt,
 *				 REAL rminfrac, short *pnc)
 *
 *	Arguments:	(x[], y[]) = Used to hold the coordinates of a random
 *				set of points (generated by this routine)
 *			(wid, hgt) = Dimensions of box in which the generated
 *				random points will be constrained to lie
 *			rminfrac = Parameter to determine degree of randomness
 *				of distribution of points (between 0 and 1)
 *			*pnc	= requested number of random points
 *
 *	Return value:	none
 *
 *	Action:		Generates a random 2-d distribution of points within
 *			a box of size [-wid/2,wid/2]x[-hgt/2,hgt/2].  The
 *			degree of randomness can be tuned using the parameter
 *			`rminfrac'.  A value of 0.0 gives complete
 *			randomness and a value of 1.0 corresponds to perfect
 *			(hexagonal) ordering of points.  The variable
 *			`rminfrac' expresses the minimum distance between
 *			particles as a fraction of the value for hexagonal
 *			close packing.
 *
 *			The number of points requested is passed in `*pnc'
 *			and the generated points are put into the arrays
 *			(x[], y[]).  Note that it is impractical to request
 *			distributions with order `rminfrac' much greater
 *			than about 0.6
 *
 *	Acknowledgements: Adapted from a routine by J. P. Kermode.
 *
 *****************************************************************************/
void vorspray(x, y, wid, hgt, rminfrac, pnc)
REAL *x, *y, wid, hgt, rminfrac;
short *pnc;
{
  REAL d, rmin, x1, y1, x2, y2, ran3(), linlen();
  int rseed;
  short i, j, nt=0, np=0;
  boolean flag;
  void plerror(), vortrans();
	struct timeval now;
  rseed= (int) svorseed;
  /* 'rminfrac' expresses the min radius as fraction of the seperation */
  /* required for hexagonal close packing of circles. */
  rmin=rminfrac*sqrt(8.0*wid*hgt/((*pnc)*sqrt(27.0)));
  if (*pnc<1) { plerror("zero cells requested in routine vorspray"); return; }
  if (rmin>wid || rmin>hgt) {
    *pnc=0; plerror("hard disc too large in routine vorspray");
    return;
  }
  if (rseed!=0) {
#ifndef FIXEDSEED
		gettimeofday(&now,NULL);
		rseed = now.tv_usec;
#endif
		ran3(-rseed);
		rseed=0; svorseed=0.0; }
  x[0]=(ran3(0)-0.5)*wid; y[0]=(ran3(0)-0.5)*hgt; np++;
  while (np<*pnc) {
    nt=0;
    do {
      flag=TRUE;
      x1=(ran3(0)-0.5)*wid; y1=(ran3(0)-0.5)*hgt; nt++;
      for (i=0; i<np; i++) {
        vortrans(x1,y1,x[i],y[i],&x2,&y2);
        if ((d=linlen(x2,y2,x[i],y[i])) < rmin) { flag=FALSE; break; }
      }
    } while (!flag && nt<10000);
    if (nt>=10000) {
      *pnc=np;
      plerror("too many attempts to place a point in routine vorspray");
      return;
    }
    x[np]=x1; y[np]=y1; np++;
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *    V O R D E L A U N E Y ( )     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	boolean vordelauney()
 *
 *	Arguments:	none
 *
 *	Return value:	Returns FALSE if the construction gets stuck
 *
 *	Action:		This routine carries out a Delauney triangulation
 *			upon the random set of points in the periodic box (as
 *			obtained from the routine `vorspray()' ).  A Delauney
 *			triangulation consists of joining each point to its
 *			nearest neighbours in such a way that the plane is
 *			completely covered by a perfect tiling of
 *			*non-overlapping* triangles.  This tiling of triangles
 *			is then the dual graph of a dry froth.  If we carry
 *			out the `duality' transformation from a network of
 *			triangles to a network of cells then we get the
 *			following mapping:
 *				triangle vertex <---> cell centre
 *						       = (cx[i], cy[i])
 *				triangle centre <---> cell vertex
 *					               = (vx[i], vy[i])
 *				no. of vertex nbrs. <---> no. of sides of cell
 *				no. of triangle vertices <---> no. of cells
 *								= `nc'
 *				no. of triangles <---> no. of cell vertices
 *							= `nv'
 *			Note that the names given to the variables we use
 *			reflect the structure of the dual (dry froth) network.
 *
 *			THE `TRIANG[][][]' ARRAY:
 *			The array used to hold the topological structure of
 *			the triangulation is `triang[][][]' (this does not
 *			appear explicitly here but is accessed via the two
 *			subroutines `dltri1()' and `dltri2()' ) and
 *			geometrical information is stored in the arrays
 *			(vx[], vy[]) which store the centre points of the
 *			triangles.  The format of the array `triang[][][]' is
 *			as follows:
 *				triang[triangle_index][nbr_list][info_type]
 *			where the `triangle_index' (in the range 0..nv)
 *			selects a particular triangle, `nbr_list' (in the
 *			range `k = 0, 1, 2') selects one of its three
 *			neighbours and `info_type' is one of:
 *				0  => index of `nbr_list'th  nbr.
 *				1  => periodic index of `nbr_list'th  nbr.
 *			(the business of periodic indices is a bit tricky,
 *			essentially telling you how to put the nbr. into the
 *			same periodic box as the first triangle -- consult the
 *			suite of routines `vorvnbrxy()', `PERFN()', `PERX()'
 *			and `PERY()' for more details).  Thus, in summary,
 *			array `triang' essentially tells you who your
 *			neighbours are.
 *
 *			PERIODIC INDICES:
 *			Throughout the Delauney routines it is important to
 *			keep in mind that an (index, periodic index) go
 *			together as a pair.  Both indices, e.g. (i, i1), are
 *			essentially needed to uniquely identify a point
 *			(where `i1' essentially tells you which periodic box
 *			the point is in).
 *
 *			THE ALGORITHM:
 *			The tricky aspect of the triangulation algorithm is
 *			that the triangles have to be chosen so as to be
 *			non-overlapping.  There is however a criterion you
 *			can use when choosing a triangle that guarantees your
 *			tiling will not overlap.  This careful `choice' of a
 *			triangle is farmed out to the two subroutines
 *			`dltri1()' and `dltri2()' (and see those routines for
 *			further details).
 *
 *			We assume the existence of these two engines
 *			`dltri1()' and `dltri2()' which perform the following
 *			tasks:
 *
 *			dltri1(&i,&i1,&j,&j1,&k,&k1) -- search for any three
 *			points, not fully enclosed within the (possibly) already
 *			partially formed network, which would make a
 *			satisfactory triangle for the tiling.  The points are
 *			returned as three (index, periodic index) pairs viz.
 *			(i, i1), (j, j1), (k, k1).
 *
 *			dltri2(i,i1,j,j1,&k,&k1,jj,jj1) -- given two points
 *			of a Delauney triangle, (i, i1) and (j, j1), find a
 *			third point (k, k1) *but* don't return with point
 *			(jj, jj1) because we already have that point.
 *
 *			The algorithm goes as follows:
 *			(i) Find a seed triangle using `dltri1()' to get three
 *			points (i,i1), (j,j1), (k,k1).
 *
 *			(ii) Construct a `wheel' using the array `wheel[][]'.
 *			Essentially we focus on point `i' (as found in (i) )
 *			and use it as a `pivot' point.  A sequence of
 *			triangles are found (using `dltri2()' ) adding spokes
 *			to the wheel.  All the points around the rim of the
 *			wheel are stored in
 *			(index= wheel[rim_index][0],
 *		 	              periodic index= wheel[rim_index][1] )
 *			with the `pivot' point `i' remaining the centre of
 *			the wheel.  This process ends when the first point
 *			of the rim (stored as (je,j1) ) equals the last
 *			point found.
 *
 *			(iii) Once a wheel-ful of triangles has been found,
 *			try using one of the points on the rim of the wheel
 *			as a new pivot point (but check it has not been
 *			used as a pivot before by looking up the array
 *			`waspivot[]').  If a pivot point is found then do
 *			step (iv) else step (v).
 *
 *			(iv) Since this next wheel overlaps with the last one,
 *			you only need to construct a partial wheel.  Mark the
 *			expected endpoint in the variables (je, je1) and then
 *			go back to the beginning of the do loop...
 *
 *			(v) Because of the periodic boundary conditions, it is
 *			perfectly possible that the untiled part can break
 *			into disconnected regions.  Then the `wheel'
 *			construction will ultimately box itself into a corner.
 *			When that happens, a pivot point is not found in step
 *			(iii) and we have to search for a new seed triangle
 *			using `dltri1()'...then go back to the beginning
 *			of the do loop...
 *
 *	Acknowledgements: Adapted from a routine by J. P. Kermode.
 *
 *****************************************************************************/
boolean vordelauney()
{
  short wheel[MWHEEL][2];
  short i, i1, j, j1, k, k1, jj, jj1, je, je1, l, ll, nwh, npv;
  boolean flag, dltri1(), dltri2(), inbox();
  void plerror();
  for (i=0; i<nc; i++) waspivot[i]=FALSE;
  npv=0;		/* ...number of pivot points used so far */
  i=i1=0;
  /*
   * Find a seed triangle to start things off...
   */
  if (!dltri1(&i,&i1,&j,&j1,&k,&k1)) {
    plerror("no seed found in routine vordelauney");
    return FALSE;
  }
  je=j; je1=j1; nwh=0;
  do {
    /*
     * Construct a `wheel' of triangles, storing the points
     * along the rim in
     * (index= wheel[rim_index][0], periodic index= wheel[rim_index][1] )
     */
    do {
      wheel[nwh][0]=j; wheel[nwh][1]=j1; nwh++;
      jj=j; jj1=j1;
      j=k; j1=k1;
      dltri2(i,i1,j,j1,&k,&k1,jj,jj1);
      if (nv>=(nc+nc)) return TRUE;
    } while ((k>=0) && (k!=je || k1!=je1));
    waspivot[i]=TRUE; npv++;
    j=i; j1=i1; flag=FALSE;
    /*
     * Search the rim of the wheel for a new pivot point...
     */
    for (l=1; (l<nwh-1) && !flag; l++) {
      i=wheel[l][0]; i1=wheel[l][1];
      if (!waspivot[i] && inbox(i,i1)) {
        ll=(l+1)%nwh;
        k=wheel[ll][0]; k1=wheel[ll][1];
        ll=(l+nwh-1)%nwh;
        je=wheel[ll][0]; je1=wheel[ll][1];
        nwh=0;
        wheel[nwh][0]=je; wheel[nwh][1]=je1; nwh++;
        flag=TRUE;
        break;
      }
    }
    /*
     * If a new pivot point could not be found on the rim
     * then cast for a completely new seed triangle to
     * get things going again...
     */
    if (!flag) {
      if (!dltri1(&i,&i1,&j,&j1,&k,&k1)) {
        plerror("no seed found in routine vordelauney");
        return FALSE;
      }
      if (nv>=(nc+nc)) return TRUE;
      je=j; je1=j1; nwh=0;
    }
  } while (npv<nc);
  return FALSE;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  D L T R I 1 ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	boolean dltri1(short *i, *i1, *j, *j1, *k, *k1)
 *
 *	Arguments:	(*i, *i1), (*j, *j1), (*k, *k1) = Three pairs of
 *			indices in the format (index, periodic index) are
 *			returned by this routine.  These indices identify
 *			three points which form a Delauney triangle.
 *
 *	Return value:	TRUE if a triangle is found, otherwise FALSE.
 *
 *	Action:		Used to find the three coordinates of a `seed'
 *			Delauney triangle -- either to start the network off
 *			or to get it going again if it gets stuck.  It returns
 *			an (index, periodic index) pair for the three apices
 *			of a triangle satisfying the Delauney condition.
 *			The condition states that when we draw a circumcircle
 *			around the triangle there should be *no* other points
 *			lying inside this circle.
 *
 *			There are two different searches attempted.  First of
 *			all a search is made of the existing network to find
 *			a triangle containing a point (*i, *i1) which has not
 *			yet been used as a pivot point (which implies that it
 *			lies on the boundary of the triangle network).  If
 *			this succeeds then it is by far the quicker method.
 *
 *			If the first method fails or there is no existing
 *			network, then a Delauney triangle has to be found from
 *			scratch.  This involves a search using a hefty number
 *			of nested loops! (see accompanying comments).
 *
 *			This routine guarantees that `waspivot[*i]=FALSE'
 *			whenever a triangle is returned.
 *
 *	Acknowledgements: Adapted from a routine by J. P. Kermode.
 *
 *****************************************************************************/
boolean dltri1(i, i1, j, j1, k, k1)
short *i, *i1, *j, *j1, *k, *k1;
{
  boolean flag, tricen();
  short t, l, ll, l1, jo[9], ko[9], i3, j3, k3, m, m1, mm1, mo[9];
  REAL x1, y1, x2, y2, x3, y3, xc, yc, r1, r2, wid=boxwid, hgt=boxhgt,
    rr, xm, ym, vorrad, rrlim;
  REAL linlen2();
  void plerror(), vororder(), ptbox();
  /* vorrad= VORRAD*sqrt(boxwid*boxhgt/((REAL) nc)); */
  /*
   *  Set `vorrad' to be roughly equal to the distance between two
   *  points in the random network...
   */
  vorrad=1.1*srfrac*sqrt(8.0*boxwid*boxhgt/(nc*sqrt(27.0)));
  /*
   * Upper limit for the search area (rather large)...
   */
  rrlim=4000.0*vorrad*vorrad;
  /*
   * Lists the optimum order in which to search neighbouring periodic boxes
   * for a third point in `mo[]', when the first two points are in the
   * central periodic box...
   */
  vororder(0,0,mo);
  /*
   * First do a search of the triangle network to see if you can
   * find a triangle plus a point on this triangle which has not been used
   * before as a pivot point.  If we succeed then we're done...
   */
  if (nv!=0) {
    for (t=0; t<nv; t++) {
      for (l=0; l<3; l++) {
        *i=triang[t][l][0];
        if (!waspivot[*i]) {
          *i1=triang[t][l][1];
          ll=(l+1)%3;
          *j=triang[t][ll][0]; *j1=triang[t][ll][1];
          ll=(l+2)%3; *k=triang[t][ll][0]; *k1=triang[t][ll][1];
          return TRUE;
        }
      }
    }
    /*
     * ...else take any ol' point (*i, *i1) which has not been used as a pivot point
     * before...
     */
    for (*i=0; *i<nc; (*i)++) { if (!waspivot[*i]) break; }
  }
  *i1=0;
  rr=vorrad*vorrad;	/* ...initial search radius */
  x1=cx[*i]; y1=cy[*i];
  /*
   * Put optimised order to search periodic boxes into `jo[]' given that
   * the first point has periodic index *i1
   */
  vororder(*i1,*i1,jo);
  /*
   * `do' loop for ever increasing search area `rr'
   * until a Delauney triangle is found...
   */
  do {
    /*
     * Loop over all possible indices (*j, *j1) such that
     *   (i)	(*i, *i1) not equiv (*j, *j1),
     *   (ii)	*j has not been used as a pivot point
     *   (iii)	point *j lies within the search circle
     * then given this (*j, *j1), do more nested loops to find (*k, *k1) s.t.
     *   (i)	(*k, *k1) distinct from both (*i, *i1) and (*j, *j1),
     *   (ii)	*k has not been used as a pivot point
     *   (iii)	point *k also lies within the search circle
     * Now that we have the prospective triangle given by points
     * (*i, *i1), (*j, *j1), (*k, *k1) we have to test for the Delauney
     * condition:  This state that if we draw a circle through the apices
     * of the triangle then there must be no other points within this circle.
     * To test this we loop through all other points indexed by (*m, *m1) s.t.
     *   (i)	(*m, *m1) distinct from the three apices of the triangle
     * and we require that for *all* of these points that they lie outside the
     * circle of centre (xc, yc) and radius `r1', otherwise set flag=FALSE.
     *
     * If flag remains TRUE then the triangle has passed the Delauney criterion
     * and the search has been successful.
     */
    for (l=0, *j1=jo[0]; l<9; l++, *j1=jo[l]) {
      for (*j=0; *j<nc; (*j)++) {
        if ((*i!=*j || *i1!=*j1) && !waspivot[*j]) {
          x2=cx[*j]+wid*PERX(*j1); y2=cy[*j]+hgt*PERY(*j1);
          if (linlen2(x1,y1,x2,y2)<=rr) {
            vororder(*i1,*j1,ko);
            for (l1=0, *k1=ko[0]; l1<9; l1++, *k1=ko[l1]) {
              for (*k=0; *k<nc; (*k)++) {
                if ((*i!=*k || *i1!=*k1) && (*j!=*k || *j1!=*k1) &&
                  !waspivot[*k]) {
                  x3=cx[*k]+wid*PERX(*k1); y3=cy[*k]+hgt*PERY(*k1);
                  if (linlen2(x1,y1,x3,y3)<=rr) {
                    if (tricen(x1,y1,x2,y2,x3,y3,&xc,&yc)) {
                      r1=linlen2(xc,yc,x1,y1);
                      flag=TRUE;
                      for (m=0; (m<nc) && flag; m++)
                        for (m1=mo[mm1=0]; (mm1<9) && flag; m1=mo[++mm1])
                          if ((*i!=m || *i1!=m1) && (*j!=m || *j1!=m1) &&
                            (*k!=m || *k1!=m1)) {
                            xm=cx[m]+wid*PERX(m1); ym=cy[m]+hgt*PERY(m1);
                            r2=linlen2(xc,yc,xm,ym);
                            if (r2<r1) flag=FALSE;
                          }
                      if (flag) {
			/*
			 * Delauney triangle has been successfully found so
			 * add it to the network here...
			 */
			/*
			 * A slightly tricky operation with periodic indices:
			 * If coordinates (xc, yc) are not in the central periodic
			 * box then `ptbox()' translates them into it, *and* at
			 * the same updates i3, j3, k3 so that they will translate
			 * the respective points into the vicinity of the new values
			 * of the coordinates (xc, yc).
			 */
                        i3= *i1; j3= *j1; k3= *k1; ptbox(&xc,&yc,&i3,&j3,&k3);
			/*
			 * Record the centre point of the triangle.  Note that ultimately
			 * these coordinates will be the vertices of the dry foam network
			 */
                        vx[nv]=xc; vy[nv]=yc;
                        triang[nv][0][0]= *i; triang[nv][0][1]=i3;
                        triang[nv][1][0]= *j; triang[nv][1][1]=j3;
                        triang[nv][2][0]= *k; triang[nv][2][1]=k3;
                        nv++;
                        return TRUE;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    /*
     * Crank up the search radius and try again...
     */
    rr *= 1.1;
  } while (rr<rrlim);
  plerror("failure of routine dltri1");
  return FALSE;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  D L T R I 2 ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	boolean dltri2(short i, i1, j, j1, *k, *k1, jj, jj1)
 *
 *	Arguments:	(i, i1), (j, j1)  = the (index, periodic index) pairs
 *			for the first two apices of a Delauney triangle.
 *
 *			(*k, *k1)	= the (index, periodic index)
 *			pair for the third apex of a Delauney triangle.
 *
 *			(jj, jj1)	= an (index, periodic index) pair
 *			which is not to be returned in the result (*k, *k1)
 *			because this third apex has already been found.
 *
 *	Return value:	TRUE if successful, otherwise returns FALSE.
 *
 *	Action:		This routine is given two adjacent points indexed by
 * 			(i, i1) and (j, j1) on the triangulation network.  It
 * 			then finds a third point (*k, *k1) such that all three
 * 			points give a Delauney triangle, satisfying the
 * 			requisite conditions.  The extra indices (jj, jj1)
 * 			specify a point which you do *not* want returned in
 * 			(*k, *k1) because you already have this point (in a
 * 			tiling of triangles, a given edge will be adjacent to
 * 			exactly two triangles).
 *
 * 			Two different attempts can be made to find a suitable
 * 			triangle.  Firstly, the existing network of triangles
 * 			is searched to see if a suitable triangle can be
 * 			found there, with third vertex *k not equal to jj.
 *
 * 			Otherwise, a new triangle must be found.  A search is
 * 			made for a point (*k, *k1) which completes a Delauney
 * 			triangle with (i, i1) and (j, j1).
 *
 *	Acknowledgements: Adapted from a routine by J. P. Kermode.
 *
 *****************************************************************************/
boolean dltri2(i, i1, j, j1, k, k1, jj, jj1)
short i, i1, j, j1, *k, *k1, jj, jj1;
{
  short ix, iy, ix1, iy1, l, m, n, t, ko[9], m1, mm1, mo[9], i3, j3, k3;
  REAL rr, r1, r2, wid, hgt, x1, y1, x2, y2, x3, y3, xc, yc, xm, ym, vorrad,
       rrlim;
  REAL linlen2();
  boolean flag, tricen();
  void vororder(), ptbox(), plerror();
  /* vorrad= VORRAD*sqrt(boxwid*boxhgt/((REAL) nc)); */
  /*
   *  Set `vorrad' to be roughly equal to the distance between two
   *  points in the random network...
   */
  vorrad=1.1*srfrac*sqrt(8.0*boxwid*boxhgt/(nc*sqrt(27.0)));
  /*
   * Upper limit for the search area (rather large)...
   */
  rrlim=4000.0*vorrad*vorrad;
  /*
   * Lists the optimum order in which to search neighbouring periodic boxes
   * for a third point in `mo[]', when the first two points are in the
   * central periodic box...
   */
  vororder(0,0,mo);
  /*
   * First of all try and find a triangle which lies within the existing
   * Delauney network (if the network already partly exists)...
   */
  if (nv>=1) {
    /*
     * Note that `PERX()' and `PERY()' are macros which turn a periodic
     * index into the X and Y component of an integer vector.  This
     * integer vector tells you how many periodic boxes to move the point
     * up/down or left/right.
     *
     * Here we are recording the periodic box difference-vector, pointing
     * from `i' to `j', to use for comparisons later on...
     */
    ix=PERX(j1)-PERX(i1); iy=PERY(j1)-PERY(i1);
    /*
     * Loop over all triangles `t', and neighbour points on that triangle
     * given by `l = 0, 1, 2' and `m = 0, 1, 2' -- see if this triangle
     * fits (i, i1), (j, j1) and (*k, *k1)...
     */
    for (t=0; t<nv; t++)
      for (l=0; l<3; l++)
        for (m=0; m<3; m++)
          if (m!=l)
            if (i==triang[t][l][0] && j==triang[t][m][0])
              if ((ix==(PERX(triang[t][m][1])-PERX(triang[t][l][1]))) &&
                (iy==(PERY(triang[t][m][1])-PERY(triang[t][l][1])))) {
		/*
		 * ...at this stage `i' and`j' are identical to two of the
		 * points on the triangle `t' ...
		 */
                n=3-(l+m);  /* nbr. index of the third point on the triangle */
                *k=triang[t][n][0]; *k1=triang[t][n][1];
                if (*k!=jj
                   || (PERX(*k1)-PERX(jj1))!=(PERX(triang[t][l][1])-PERX(i1))
                   || (PERY(*k1)-PERY(jj1))!=(PERY(triang[t][l][1])-PERY(i1))) {
		   /*
		    * If point (*k, *k1) is not identical to (jj, jj1) then we
		    * have successfully found a triangle.
		    *
		    * Note that we cannot compare *k1 and jj1 directly, since we
		    * are not sure in which periodic box the points given in the
		    * argument list are centred.  We get around this by
		    * translating the periodic indices into vectors (using
		    * PERX() and PERY() ), and then only comparing *relative*
		    * periodic vectors...
		    */
		  /*
		   * The integer vector (ix1, iy1) tells us where the periodic
		   * box of the triangle `t' lies with respect to its image
		   * given by (i, i1), (j, jj1), (*k, *k1).  (This complication
		   * comes about because we cannot assume that the given
		   * triangle lies in the same periodic box as the triangle `t'
		   * that we find). Note that (triang[t][l][0], triang[t][l][1])
		   * is the equivalent image of (i, i1).
		   */
                  ix1=PERX(triang[t][l][1])-PERX(i1);
                  iy1=PERY(triang[t][l][1])-PERY(i1);
		  /*
		   * The vector (ix1, iy1) is then used to translate *k1 to a
		   * location which is consistent with (i, i1) and (j, j1).
		   *
		   * Note that PERFN() is used to translate a periodic vector
		   * back into a periodic index...
		   */
                  *k1=PERFN(PERX(*k1)-ix1,PERY(*k1)-iy1);
                  return TRUE;
                }
              }
  }
  /* The initial search radius is `rr'... */
  rr=vorrad*vorrad; wid=boxwid; hgt=boxhgt;
  x1=cx[i]+wid*PERX(i1); y1=cy[i]+hgt*PERY(i1);
  x2=cx[j]+wid*PERX(j1); y2=cy[j]+hgt*PERY(j1);
  /*
   * Lists the optimum order in which to search neighbouring periodic boxes
   * for a third point into `ko[]', when the first two points are in the
   * boxes determined by i1 and j1.
   */
  vororder(i1,j1,ko);
  /*
   * ...if the first search failed...
   * `do' loop over steadily increasing search radii...
   */
  do {
    /*
     * Loop over all indices (*k, *k1)...
     */
    for (l=0, *k1=ko[l]; l<9; l++, *k1=ko[l])
      for (*k=0; *k<nc; (*k)++)
        if ((i!=*k || i1!=*k1) && (j!=*k || j1!=*k1) && (jj!=*k || jj1!=*k1)
                && !waspivot[*k]) {
	  /*
	   * If (*k, *k1) is distinct and has not yet been used as pivot...
	   */
          x3=cx[*k]+wid*PERX(*k1); y3=cy[*k]+hgt*PERY(*k1);
          if (linlen2(x1,y1,x3,y3)<=rr) {
    /*
     * ...and lies within search area, then do the test for the Delauney
     * condition:  This states that if we draw a circle through the apices
     * of the triangle then there must be no other points within this circle.
     * To test this we loop through all other points indexed by (*m, *m1) s.t.
     *   (i)	(*m, *m1) distinct from the three apices of the triangle
     * and we require that for *all* of these points that they lie outside the
     * circle of centre (xc, yc) and radius `r1', otherwise set flag=FALSE.
     */
            flag=FALSE;
            if (tricen(x1,y1,x2,y2,x3,y3,&xc,&yc)) {
              r1=linlen2(xc,yc,x1,y1);
              flag=TRUE;
              for (m1=mo[mm1=0]; (mm1<9) && flag; m1=mo[++mm1])
                for (m=0; (m<nc) && flag; m++)
                  if ((i!=m || i1!=m1) && (j!=m || j1!=m1) && (*k!=m || *k1!=m1)) {
                    xm=cx[m]+wid*PERX(m1); ym=cy[m]+hgt*PERY(m1);
                    r2=linlen2(xc,yc,xm,ym);
                    if (r2<r1) flag=FALSE;
                  }
            }
            if (flag) {
	      /*
	       * Delauney triangle has been successfully found so
	       * add it to the network here...
	       */
	      /*
	       * A slightly tricky operation with periodic indices:
	       * If coordinates (xc, yc) are not in the central periodic
	       * box then `ptbox()' translates them into it, *and* at
	       * the same updates i3, j3, k3 so that they will translate
	       * the respective points into the vicinity of the new values
	       * of the coordinates (xc, yc).
	       */
              i3=i1; j3=j1; k3= *k1; ptbox(&xc,&yc,&i3,&j3,&k3);
	      /*
	       * Record the centre point of the triangle.  Note that ultimately
	       * these coordinates will be the vertices of the dry foam network
	       */
              vx[nv]=xc; vy[nv]=yc;
              triang[nv][0][0]=i; triang[nv][0][1]=i3;
              triang[nv][1][0]=j; triang[nv][1][1]=j3;
              triang[nv][2][0]= *k; triang[nv][2][1]=k3;
              nv++;
              return TRUE;
            }
          }
        }
    /*
     * Crank up the search radius and try again...
     */
    rr *= 1.1;
  } while (rr<rrlim);
  *k = -1;  /* ...value indicates error return */
  plerror("failure of routine dltri2");
  return FALSE;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *    V O R D R Y T O P O L ( )     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	vordrytopol()
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		After the Delauney triangulation has been found it is
 *			necessary to convert the triangle topology, stored in
 *			`triang[][][]', into the topology of its dual graph
 *			to get a dry froth network.  Under the duality
 *			transformation we have some equivalences:
 *
 *			centre of triangle (vx[], vy[]) <--> new vertex coord
 *			index of triangle `i' <--> index of new vertex
 *			apex of triangle  <-->  centre of new cell
 *			index of triangle apex <--> index of new cell
 *
 *			This routine has essentially to initialise the
 *			matrices `vorvnbr[][]', `vorcadj[][]' and
 *			`vorvper[][]' which define the topology of the dry
 *			froth network.  The main loop of this routine runs
 *			through all possible pairs of triangles `i' and `j'
 *			in an attempt to find a pair which share one edge in
 *			common (which means we also have to loop over all of
 *			the edges `k = 0, 1, 2' of each of the triangles).
 *			Once an adjacent pair of triangles has been found we
 *			are able to fill in some information about the topology
 *			of the dry froth.
 *
 *	Acknowledgements: Adapted from a routine by J. P. Kermode.
 *
 *****************************************************************************/
void vordrytopol()
{
  short i, j, k, k1, px, py, px1, py1, (*t)[MBORD][3][2];
  /*
   * The matrix ik[][] gives a simple way of selecting neighbour indices
   * `k = 0, 1, 2'.  If we are given the k'th index, then the *other* two
   * indices are the `ik[k][0]'th and the `ik[k][1]'th  (out of the set
   * {0, 1, 2} ).
   */
  static short ik[3][2]={ {1, 2}, {2, 0}, {0, 1} };
  REAL x1, y1, x2, y2, x3, y3, cross;
  void vorvnbrxy(), plerror();
  t= &triang;  /* Abbreviation for the `triang[][][]' array */
  for (i=0; i<nv; i++) for (k=0; k<3; k++) vorvnbr[i][k]= -1;
  /*
   * Convert the triangle topology into a dry froth topology.
   * Note that the centre points of the triangle (vx[i], vy[i])
   * will correspond to the new vertices of the dry froth.
   * Therefore to say that "the k'th neighbour of triangle `i' is
   * triangle `j' ", is equivalent to saying that the k'th neighbour
   * of new vertex `i' is new vertex `j'.
   */
  /*
   * Loop over all triangles (vertices) `i' and all triangles (vertices) `j'.
   * They each loop over their three neighbours indices `k' and `k1'
   * respectively, to see if they are each others neighbour.
   */
  for (i=0; i<nv-1; i++)
    for (k=0; k<3; k++)
      if (vorvnbr[i][k]<0)
        for (j=i+1; j<nv; j++)
          for (k1=0; k1<3; k1++)
            if (vorvnbr[j][k1]<0)
              if ((*t)[i][ik[k][0]][0]==(*t)[j][ik[k1][0]][0]) {
                if ((*t)[i][ik[k][1]][0]==(*t)[j][ik[k1][1]][0]) {
		  /*
		   * ...implies two points on the boundary of triangle `i' and
		   * two points on the boundary of triangle `j' match
		   * up (in such a way that implies `j' is the k'th nbr. of
		   * `i' and `i' is the k1'th nbr. of `j')
		   */
		  /*
		   * Calculate the periodic vector difference between two
		   * points on each of the triangles which are equivalent
		   * ==> vector (px, py) tells by how many periodic boxes
		   * the triangles (new vertices) must be shifted so as to
		   * line up.
		   *
		   * Note that calculating (px1, px2) as well is redundant
		   * (since it should be equal to (px, py) )
		   * and is only done as an elementary means of error checking.
		   */
                  px=PERX((*t)[j][ik[k1][0]][1])-PERX((*t)[i][ik[k][0]][1]);
                  py=PERY((*t)[j][ik[k1][0]][1])-PERY((*t)[i][ik[k][0]][1]);
                  px1=PERX((*t)[j][ik[k1][1]][1])-PERX((*t)[i][ik[k][1]][1]);
                  py1=PERY((*t)[j][ik[k1][1]][1])-PERY((*t)[i][ik[k][1]][1]);
                  if (px==px1 && py==py1) {
                    vorvnbr[i][k]=j; vorvnbr[j][k1]=i;
		    /*
		     * The apices of the triangles are the dual of the cells
		     * of the dry froth network...so use their indices now as
		     * cell indices.
		     */
                    vorcadj[i][k]=(*t)[i][k][0]; vorcadj[j][k1]=(*t)[j][k1][0];
                    vorvper[i][k]=PERFN(-px,-py);
                    vorvper[j][k1]=PERFN(px,py);
                    if (abs(px)>1 || abs(py)>1)
                      plerror("topological error in routine vordrytopol");
                  }
                }
              }
              else if ((*t)[i][ik[k][0]][0]==(*t)[j][ik[k1][1]][0]) {
                if ((*t)[i][ik[k][1]][0]==(*t)[j][ik[k1][0]][0]) {
		  /*
		   * ...two points on the boundary of triangle `i' and
		   * on the boundary of triangle `j' match up contrariwise.
		   * (in such a way that implies `j' is the k'th nbr. of
		   * `i' and `i' is the k1'th nbr. of `j')
		   */
		  /* ...see comments for block above ... */
                  px=PERX((*t)[j][ik[k1][1]][1])-PERX((*t)[i][ik[k][0]][1]);
                  py=PERY((*t)[j][ik[k1][1]][1])-PERY((*t)[i][ik[k][0]][1]);
                  px1=PERX((*t)[j][ik[k1][0]][1])-PERX((*t)[i][ik[k][1]][1]);
                  py1=PERY((*t)[j][ik[k1][0]][1])-PERY((*t)[i][ik[k][1]][1]);
                  if (px==px1 && py==py1) {
                    vorvnbr[i][k]=j; vorvnbr[j][k1]=i;
                    vorcadj[i][k]=(*t)[i][k][0]; vorcadj[j][k1]=(*t)[j][k1][0];
                    vorvper[i][k]=PERFN(-px,-py);
                    vorvper[j][k1]=PERFN(px,py);
                    if (abs(px)>1 || abs(py)>1)
                      plerror("topological error in routine vordrytopol");
                  }
                }
              }
  /*
   * This loop essentially enforces the convention on how neighbouring
   * vertices are indexed...
   *
   * The neighbours of vertex `i' are indexed by `k = 0, 1, 2'.
   * It is convenient to adopt a convention whereby these
   * neighbours are always referred to in anti-clockwise order
   * as you go through the index `k = 0, 1, 2'.  This is exactly what
   * this loop does:  It loops through all vertices and for each vertex
   * it constructs a triangle consisting of the three neighbouring
   * vertex coordinates.  A cross product is used to establish the sense
   * of this triangle.  If this cross product is opposite to the desired
   * sense, then two of the vertices are swapped over (achieved by
   * swapping the information in the topological arrays).
   *
   * Note that here, for the first time, we use the routine
   * `vorvnbrxy(i, k, &x, &y)'.  This is a useful routine which
   * returns the coordinates (x, y) of the `k'th neighbour of vertex `i'.
   * This humble task is just a bit messier than you might expect on account
   * of the periodic boundary conditions!  (you want the neighbours all in
   * same periodic box after all).  Much heartache is avoided by using
   * `vorvnbrxy()' as a black box...
   */
  for (i=0; i<nv; i++) {
    for (k=0; k<3; k++)
      if (vorvnbr[i][k]<0)
        plerror("topology incomplete in routine vordrytopol");
      else {
        vorvnbrxy(i,0,&x1,&y1);
        vorvnbrxy(i,1,&x2,&y2);
        vorvnbrxy(i,2,&x3,&y3);
        x2 -= x1; y2 -= y1;
        x3 -= x1; y3 -= y1;
        cross=x2*y3-x3*y2;
        if (cross<0.0) {
          j=vorvnbr[i][1]; vorvnbr[i][1]=vorvnbr[i][2]; vorvnbr[i][2]=j;
          j=vorcadj[i][1]; vorcadj[i][1]=vorcadj[i][2]; vorcadj[i][2]=j;
          j=vorvper[i][1]; vorvper[i][1]=vorvper[i][2]; vorvper[i][2]=j;
        }
      }
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    * V O R D E L T A S E A R C H ( )     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL vordeltasearch()
 *
 *	Arguments:	none
 *
 *	Return value:	An appropriate size `delta' for the Plateau borders,
 *			which needs to be known in order to carry out the
 *			conversion of a dry froth to one with Plateau borders
 *
 *	Action:		This routine operates on a dry foam network and
 *			searches for an appropriate size of Plateau border
 *			`delta' which can easily be superimposed upon the
 *			existing vertices.  When the 3-sided borders are
 *			introduced on top of the vertices, their three apices
 *			will be moved back from the idealised vertex by a
 *			distance `delta'.  A problem arises if two of the
 *			original vertices are so close together that the
 *			newly introduced Plateau borders would overlap.  This
 *			routine is essentially intended to cope with this
 *			difficulty.
 *
 *			The strategy is to set `delta' initially to the value
 *			`smaxbfrac * f' where `f' is roughly a typical cell
 *			diameter.  If this value of `delta' brings two Plateau
 *			borders too close together than `delta' is reduced
 *			just enough to avoid the clash.  However, since it is
 *			difficult to equilibrate very small Plateau borders
 *			we do not allow `delta' to become arbitrarily small.
 *			There is a lower limit of `sminbfrac * f'.  If two
 *			vertices are found to be closer than allowed by this
 *			limit than they are nudged apart by this subroutine
 *			to get them out of danger.
 *
 *****************************************************************************/
REAL vordeltasearch()
{
  short i, j, j1, j2, k, k1, k2;
  REAL d, d1, d2, dmin, dmax, dx, dy, dx1, dy1, dx2, dy2, xy, x1, y1,
       x10, y10, x11, y11, x12, y12, f, tn, dxx, dyy, dd;
  REAL fsign(), sqf();
  void vorvnbrxy(), plerror();
  f=sqrt(boxwid*boxhgt/((REAL) nc));
  dmax=smaxbfrac*f; dmin=sminbfrac*f;
  d=dmax;
  for (i=0; i<nv; i++) {
    for (k=0; k<3; k++) {
      if (vorvnbr[i][k]>i) {
        x1=vx[i]; y1=vy[i];
        vorvnbrxy(i,k,&x10,&y10);
        dx=x10-x1; dy=y10-y1;
        if (fabs(dx)<4.0*d && fabs(dy)<4.0*d) {
          if ((xy=0.25*(dx*dx+dy*dy))<d*d) {
            d=sqrt(xy);
            if (d<dmin) {
              f=dmin/d-1.0;
              j=vorvnbr[i][k];
              vx[i] -= 0.5*f*dx+0.2*dmin*fsign(dx);
              vy[i] -= 0.5*f*dy+0.2*dmin*fsign(dy);
              vx[j] += 0.5*f*dx+0.2*dmin*fsign(dx);
              vy[j] += 0.5*f*dy+0.2*dmin*fsign(dy);
              d=dmin;
              plerror("delta too small in routine vdeltasearch");
              /* now check that it is not now too close to its nbrs */
              k1=(k+1)%3; k2=(k+2)%3;
              j1=vorvnbr[i][k1]; j2=vorvnbr[i][k2];
              vorvnbrxy(i,k1,&x11,&y11); dx1=x11-x1; dy1=y11-y1;
              vorvnbrxy(i,k2,&x12,&y12); dx2=x12-x1; dy2=y12-y1;
              if ((d1=0.5*sqrt(dx1*dx1+dy1*dy1)) < dmin) {
                /* the new co-ords for point 'i' are got as the apex of an */
                /* isosceles triangle with equal sides of length 2.4*dmin */
                dxx=x11-x10; dyy=y11-y10; dd=0.5*sqrt(dxx*dxx+dyy*dyy);
                tn=sqrt(sqf(2.4*dmin/dd)-1.0);
                tn *= fsign(dx*dy1-dx1*dy);
                vx[i]=x10+0.5*(dxx-tn*dyy);
                vy[i]=y10+0.5*(dyy+tn*dxx);
              }
              else if ((d2=0.5*sqrt(dx2*dx2+dy2*dy2)) < dmin) {
                dxx=x10-x12; dyy=y10-y12; dd=0.5*sqrt(dxx*dxx+dyy*dyy);
                tn=sqrt(sqf(2.4*dmin/dd)-1.0);
                tn *= fsign(dx2*dy-dy2*dx);
                vx[i]=x12+0.5*(dxx-tn*dyy);
                vy[i]=y12+0.5*(dyy+tn*dxx);
              }
            }
          }
        }
      }
    }
  }
  return d;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *V O R P L A T E A U T O P O L ( ) * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	vorplateautopol(REAL d)
 *
 *	Arguments:	d	= distance from the dry froth vertex to the
 *				new surrounding trio of vertices which
 *				define the 3-sided Plateau border
 *
 *	Return value:	none
 *
 *	Action:		This subroutine converts a dry foam network into a
 *			Plateau border foam network.  It works by decorating
 *			each vertex of the dry foam network with a *small*
 *			three-sided Plateau border.  The initial size of the
 *			Plateau borders is determined by the argument `d'.
 *
 *			This routine first creates an entirely new topology
 *			for the network.  The variables associated with the
 *			dry topology (`vorvnbr[][]', `vorvper[][]' and
 *			`vorcadj[][]') are replaced by the new variables
 *			`vnbr', `vper' and `cadj' which describe the Plateau
 *			network.  In the new network there are three times as
 *			many vertices and these new vertices are indexed by
 *			the local variable `newi'.  The old, dry vertices are
 *			indexed by `i'.  If we fix our attention on old
 *			vertex `i' for a moment -- it has three old neighbours
 *			which we may index by `k = 0, 1, 2' (a standard
 *			convention is adopted whereby the three neighbours
 *			surrounding the old vertex are indexed by `k' in
 *			*anti-clockwise* order always).
 *
 *			Similarly, the new vertices also happen to have three
 *			neighbours each and these can be indexed by an index
 *			we will refer to as `newk' (though it is not mentioned
 *			by name in this subroutine).  You can easily convince
 *			yourself by sketching the topology (or referring to
 *			the documentation) that one of these neighbours is
 *			always connected along a cell-cell arc while the other
 *			two neighbours are connected via border-cell arcs.
 *			(A `cell-cell' arc is an arc which has a cell on
 *			either side of it, while a `border-cell' arc has a
 *			Plateau border on one side and a cell on the other).
 *			In order to remember which arc is which, the
 *			neighbours are indexed in a standard order.  The
 *			neighbour which lies at the other end of the cell-cell
 *			arc is indexed by `newk = 0', and then in
 *			anti-clockwise order the next two neighbours are
 *			indexed by `newk = 1' and `newk = 2'.  The indices
 *			`newi' and `newk' are used by the arrays `vnbr' and
 *			`vper', so we can summarise their use:
 *
 *				vnbr[newi][newk] := index of the newk'th
 *				neighbour of vertex `newi', where nbr newk=0
 *				lies along a cell-cell arc.
 *
 *				vper[newi][newk] := info vis a vis periodic
 *				boundary conditions for nbr `newk' of
 *				vertex `newi'.
 *
 *			In the new topology the cells and Plateau borders are
 *			numbered separately.  The numbering of cells is simply
 *			the same as the old dry topology.  The Plateau borders
 *			are a new feature.  Since a Plateau border appears
 *			precisely where an old vertex used to be, the old
 *			vertex indices are used to number the new Plateau
 *			borders.  This brings us to the meaning of the array
 *			`cadj[][]'. The array `cadj[newi][cellnbr]' is used to
 *			record the indices of the Plateau border and the two
 *			cells which lie adjacent to the vertex `newi'.  A
 *			convention is adopted whereby `cellnbr = 0' referes to
 *			the Plateau border and `cellnbr = 1, 2' refer to the
 *			two subsequent cells as you turn *anti-clockwise*
 *			around the vertex `newi'.
 *
 *			Once the topology has been defined, the next loop is
 *			needed to actually place the new vertices in the
 *			correct position.  For this, the old arrays `vx[]' and
 *			`vy[]' are recycled and used with the new topology.
 *			This loop is a bit tricky (sorry).  By starting at the
 *			high index `nb' of the old coordinates and working
 *			downwards (and also working downwards from index
 *			`3*nb' of the old coordinates) it manages to avoid
 *			overwriting the old array values until they have been
 *			dealt with.  This works because of the simple linear
 *			relation between indices:  newi = 3*oldi + k .
 *
 *****************************************************************************/
void vorplateautopol(d)
REAL d;
{
  short i, k, j, j1, k1, newi, vorkindex();
  REAL dx, dy, f;
  void vnbrxy(), vorvnbrxy(), putinbox();
  nb=nv; nv *= 3;
  for (i=0; i<nb; i++) nbsides[i]=3;
  for (i=0; i<nb; i++)
    for (k=0; k<3; k++) {
      newi=3*i+k;
      vper[newi][0]=vorvper[i][k];
      vper[newi][1]=vper[newi][2]=0;
      cadj[newi][0]=i;
      cadj[newi][1]=vorcadj[i][(k+1)%3];
      cadj[newi][2]=vorcadj[i][(k+2)%3];
      j=vorvnbr[i][k]; j1=vorvper[i][k];
      k1=vorkindex(i,j,j1);
      vnbr[newi][0]=3*j+k1;
      vnbr[newi][1]=3*i+(k+1)%3;
      vnbr[newi][2]=3*i+(k+2)%3;
    }
  for (i=nb-1; i>=0; i--)
    for (k=2; k>=0; k--) {
      newi=3*i+k;
      j=vorvnbr[i][k];
      if (j>newi) vnbrxy(newi,0,&dx,&dy);
      else vorvnbrxy(i,k,&dx,&dy);
      dx -= vx[i]; dy -= vy[i];
      f=d/sqrt(dx*dx+dy*dy);
      vx[newi]=vx[i]+f*dx; vy[newi]=vy[i]+f*dy;
      putinbox(newi);
    }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  V O R O N O I ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	voronoi(short ncells)
 *
 *	Arguments:	ncells	= number of cells requested for new network
 *
 *	Return value:	none
 *
 *	Action:		This subroutine drives the other Voronoi routines.
 *			Upon finishing it should have generated an initial
 *			(unequilibrated) wet froth network with `ncell'
 *			cells.  One other basic parameter which determines
 *			the network is `rfrac' which controls the degree of
 *			disorder (see `vorspray()' for more details of this).
 *
 *			This routine is fairly self-explanatory (the gory
 *			details being concealed in the routines called
 *			from here).
 *
 *****************************************************************************/
void voronoi(ncells)
short ncells;
{
  short i, k, ii, c, nc1, csides();
  char s[80];
  REAL d, rfrac, cp1, dcp1, catot, ba;
  REAL vordeltasearch(), cellarea(), bordarea();
  void vorspray(), plerror(), vordrytopol(), vorplateautopol(), vfoamplot();
  boolean found[MBORD], arcbrk, vordelauney();
  rfrac=srfrac;
  /*
   *  Make repeated attempts (if necessary) to call first the
   *  routine `vorspray()' (which sprays `ncells' random points into the
   *  periodic box to serve as cell centres), and then the routine
   *  `vordelauney()' (which carries out a Delauney triangulation on
   *  these random cell centres).
   *
   *  The repeated attempts at these routines date from the bad old
   *  days when these routines were still buggy.  These loops are
   *  hopefully now obsolete...
   */
  do {
    do {
      for (i=0; i<1; i++) {
        nc1=ncells;
        sprintf(s,"voronoi: attempt to form with rfrac = %f",rfrac);
        plerror(s);
        vorspray(cx,cy,boxwid,boxhgt,rfrac,&nc1);
        if (nc1==ncells) break;
      }
      /* rfrac -= 0.01; */
    } while (rfrac>0 && nc1<ncells);
    if (rfrac<=0.0)
      plerror("failed to find a suitable hard disc in routine voronoi");
    nc=nc1; nv=0;
  } while (!vordelauney());
  /*
   * Now go through the various steps involved in converting the
   * network of triangles, generated by `vordelauney()' and stored in
   * `triang[][][]', to a wet foam network incorporating Plateau borders.
   *
   * vordrytopol():
   * 	triangle network `triang[][][]' ----> dry froth network
   * vordeltasearch():
   *	find a suitable size of Plateau border
   * vorplateautopol():
   *	dry froth network ----> wet froth network (Plateau borders)
   *
   * Then initialise some basic variables such as lists of valid indices,
   * `vlist[]', `clist[]', `blist[]' of vertices, cells, and borders.
   * `onv', `onc', `onb' record the `old' numbers of different objects,
   * i.e. the numbers of each object when the network was created.
   *
   * Also the cell pressures, border pressures, cell areas and numbers
   * of sides on cells are all calculated and stored in their respective
   * arrays.
   */
  vordrytopol();
  d=vordeltasearch();
  vorplateautopol(d);
  for (i=0; i<nv; i++) vlist[i]=i;
  for (i=0; i<nc; i++) clist[i]=i;
  for (i=0; i<nb; i++) blist[i]=i;
  onv=nv; onc=nc; onb=nb;
  bpav= -tan(PI/6.0)/d;
  for (i=0; i<nb; i++) bp[i]=bpav;
  for (i=0; i<nc; i++) {
    darea[i]=0.0; ncsides[i]=csides(i);
    cp[i]=0.0;
  }
  catot=0.0;
  for (i=0; i<onc; i++) found[i]=FALSE;
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    for (k=1; k<3; k++) {
      if (!found[c=cadj[i][k]]) {
        found[c]=TRUE;
        carea[c]=cellarea(i,k,&arcbrk);
        catot += carea[c];
        if (arcbrk) { plerror("failure of routine voronoi"); return; }
      }
    }
  }
  ba=(boxwid*boxhgt-catot)/nb;
  for (i=0; i<nb; i++) barea[i]=ba;
  henckyeps=0.0;
}

/********************************/
/* Miscellaneous Initialisation */
/********************************/

#define MPARAM 23

char param_tok[][20] = {"vdamp", "vvdamp", "pdamp", "bpdamp", "bprelax",
     "minbfrac", "maxbfrac", "lengthdelta", "pressuredelta", "diffuserate",
     "areasup", "bareasup", "maxdv", "maxdvv", "equilsup", "minvvlen",
     "filmwid", "cosminang", "rfrac", "notopol", "minareasup", "maxiter",
     "vorseed"};
REAL *param[] = { &svdamp, &svvdamp, &spdamp, &sbpdamp, &sbprelax, &sminbfrac,
     &smaxbfrac, &slengthdelta, &spressuredelta, &sdiffuserate, &sareasup,
     &sbareasup, &smaxdv, &smaxdvv, &sequilsup, &sminvvlen, &sfilmwid,
     &scosminang, &srfrac, &snotopol, &sminareasup, &smaxiter, &svorseed};

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *    C O N S T A N T O U T ( )     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	constantout(char tok[] )
 *
 *	Arguments:	tok[]	= a string identifying a parameter
 *
 *	Return value:	none
 *
 *	Action:		This routine prints the value of the parameter whose
 *			name is given by the string `tok[]'.  It does this
 *			with the help of the two arrays `param_tok[][]' and
 *			`*param[]' listed above.
 *
 *****************************************************************************/
void constantout(tok)
char *tok;
{
  short i;
  for (i=0; i<MPARAM; i++) {
    if (strcmp(tok,param_tok[i])==0) {
      printf("%s = %f\n",tok,*param[i]);
      return;
    }
  }
  printf("unrecognised parameter\n");
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *  S E T C O N S T A N T S ( )     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	setconstants()
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		Updates the values of a number of fundamental
 *			parameters of the froth network.  The need for this
 *			subroutine arises because as the froth evolves
 *			(possibly losing cells etc.) a number of length scales
 *			will change over time.  Since a number of features of
 *			the program are sensitive to these length scales it is
 *			necessary to update various parameters to account for
 *			this.
 *
 *			Any scale sensitive parameters exist in two forms,
 *			`name' and `sname'.  For example there is the
 *			variable `equilsup', which is used to determine
 *			convergence of equilibration, and its twin
 *			`sequilsup'.  The difference is that `sname' is the
 *			unscaled value and `name' is the scaled value.
 *			A number of the unscaled parameters `sname' may be
 *			changed via the subroutine `setsconst()'.
 *
 *			NB:  This subroutine must be called whenever any
 *			significant changes are made to the network to ensure
 *			that basic scalable parameters are kept up to date.
 *			(As the program stands, this routine is called
 *			whenever necessary).  If you need to introduce a new
 *			parameter which is somehow dependent on some length
 *			scale then you should update it *here* to ensure that
 *			it is regularly kept up to date.
 *
 *****************************************************************************/
void setconstants()
{
  REAL phifn();
  void plerror();
  if (nc==0) {
    plerror("routine setconstants called with zero cells");
    return;
  }
  cscale= sqrt(boxwid*boxhgt/((REAL) nc));
  bscale= fabs(1.0/bpav);
  volfrac=phifn();
  foamlike= (bscale/cscale < 0.1) ? TRUE : FALSE;
  /* 0.91 here denotes a typical volume fraction at which you consider the */
  /* transition from 'foam' to 'bubbles' to occur */
  lengthdelta= bscale*slengthdelta;
  cosminang=1.0-2.0*bpav*lengthdelta;
  cosminang= (cosminang<scosminang) ? cosminang : scosminang;
  /* gives cosminang sufficient clearance. */
  pressuredelta= spressuredelta/cscale;
  equilsup=bscale*sequilsup;
  filmwid= sfilmwid;
  minvvlen=bscale*sminvvlen;
  /* this dependence on border pressure is quite important to ensure */
  /* that two arcs are not so far apart that they break when co-alesced */
  diffuserate= sdiffuserate;
  areasup= min(sminareasup, sareasup*(bscale*bscale)/(cscale*cscale));
  bareasup= sbareasup*(bscale*bscale)/(cscale*cscale);
  bprelax=sbprelax*max(0.3, bscale/cscale);
/* *2.0*cscale*cscale; *40.0*((REAL) nc)/(bpav*bpav); */
  maxdv= smaxdv*cscale; /*maxdv= (maxdv<0.5*bscale) ? maxdv : 0.5*bscale;*/
  maxdvv=bscale*smaxdvv;
  notopol= (boolean) snotopol;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *S E T S C O N S T ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	setsconst(char tok[], REAL val)
 *
 *	Arguments:	tok[]	= a string naming a parameter to be set
 *			val	= the value to which the named parameter
 *				is to be set
 *
 *	Return value:	none
 *
 *	Action:		Used to set interesting parameters to a new value.
 *			This subroutine is called by the command interface
 *			to enable users to interactively change certain
 *			parameters of the program.
 *
 *			Looks up the string `tok[]' in the array
 *			`param_tok[][]' and, if found, sets the corresponding
 *			parameter (found using `*param[]') to the value `val'.
 *
 *****************************************************************************/
void setsconst(tok,val)
char *tok;
REAL val;
{
  void setconstants(), plerror();
  short i;
  for (i=0; i<MPARAM; i++) {
    if (strcmp(tok,param_tok[i])==0) { *param[i]= val; setconstants(); return; }
  }
  plerror("attempted to set an unrecognised parameter");
}
