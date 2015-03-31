#include "include.h"
/* #include "nrutil.c" */

/****************************/
/* 'Voronoi' minor routines */
/****************************/

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *V O R O R D E R ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	vororder(short i1, j1, ko[])
 *
 *	Arguments:	i1	= periodic index of first periodic box
 *			j1	= periodic index of second periodic box
 *			ko[]	= an ordered list of nine periodic boxes
 *				beginning with the periodic boxes which lie
 *				closest to `i1' and `j1'.
 *
 *	Return value:	none
 *
 *	Action:		This is icing on the cake for the Voronoi routines.
 *			It facilitate the search for a third point of a
 *			triangle in routine `dltri1()' by listing the optimum
 *			order in which to search neighbouring periodic boxes.
 *			It is just an optimisation to improve efficiency.
 *
 *			The periodic boxes which are considered in the Voronoi
 *			routine are just an array of 3x3 boxes given by
 *			integer indices (kx= {-1, 0, +1}, ky= {-1, 0, +1} ).
 *			The vectors which are encoded in the periodic indices
 *			`i1' and `j1' identify two periodic boxes out of these
 *			nine.  The purpose of this routine is to list the nine
 *			periodic boxes in the order of which are closest to
 *			the two given boxes (e.g. the first two boxes in this
 *			list will inevitably be the two boxes `i1' and `j1'
 *			themselves).  The resulting list of nine boxes are
 *			encoded as periodic indices and returned in the 
 *			array `ko[0..8]'.
 *
 *	Acknowledgements: Adapted from a routine by J. P. Kermode.
 *
 *****************************************************************************/
void vororder(i1, j1, ko)
short i1, j1, *ko;
{
  short kx, ky, k1, k2, k3, jo[9], kk[9], d[9], d1, d2;
  /*
   * The array variables have the following meanings:
   * kk[i]	= records the order of indices 0..8
   * jo[i]	= records the periodic index of the i'th box
   *		(corresponding to integer vector (kx, ky) )
   * d[i]	= contains the `distance' of the i'th periodic box (= jo[i])
   *		from the periodic box pair `i1' and `j1'.
   *
   */
  /*
   * Loop over all nine boxes (surrounding and including
   * the central box (kx= 0, ky= 0) )...
   */
  k1=0;
  for (kx= -1; kx<=1; kx++) {
    for (ky= -1; ky<=1; ky++) {
      kk[k1]=k1;
      jo[k1]=PERFN(kx,ky);
      /*
       * These next three lines effectively define the `metric'
       * for the distance from the k1'th box to the pair ( i1 and j1 )
       */
      d1=abs(PERX(i1)-kx)+abs(PERY(i1)-ky);
      d2=abs(PERX(j1)-kx)+abs(PERY(j1)-ky);
      d[k1]= (d1<d2) ? d1 : d2;
      k1++;
    }
  }
  /*
   * Do a crude bubble sort (mea culpa !) on the list of indices `kk[i]'...
   * (based on putting the distances `d[i]' into ascending order)
   */
  for (k1=8; k1>0; k1--)
    for (k2=k1-1; k2>=0; k2--)
      if (d[kk[k1]]<d[kk[k2]]) {
        k3=kk[k1]; kk[k1]=kk[k2]; kk[k2]=k3;
      }
  /*
   * Copy the sorted list into the returned array `ko[0..8]'...
   */
  for (k1=0; k1<9; k1++) ko[k1]=jo[kk[k1]];
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  T R I C E N ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	boolean tricen(REAL x1, y1, x2, y2, x3, y3, *xc, *yc)
 *
 *	Arguments:	(x1, y1)	= first apex of triangle
 *			(x2, y2)	= second apex of triangle
 *			(x3, y3)	= third apex of triangle
 *			(*xc, *yc)	= returned value of circumcentre
 *
 *	Return value:	FALSE if the three given points are *almost* exactly
 *			colinear; otherwise TRUE
 *
 *	Action:		Returns the circumcentre of the triangle given by
 *			(x1, y1), (x2, y2), (x3, y3) unless these three points
 *			are almost exactly colinear.
 *
 *****************************************************************************/
boolean tricen(x1, y1, x2, y2, x3, y3, xc, yc)
REAL x1, y1, x2, y2, x3, y3, *xc, *yc;
{
  REAL x12, y12, x23, y23, m12, m23;
  void plerror();
  x12=(x1+x2)/2.0; y12=(y1+y2)/2.0;
  x23=(x2+x3)/2.0; y23=(y2+y3)/2.0;
  m12=y2-y1; m23=y3-y2;
  if (fabs(m12)>1e-6) {
    m12 = -(x2-x1)/m12;
    if (fabs(m23)>1e-6) {
      m23= -(x3-x2)/m23;
      if (fabs(m12-m23)>fabs(1.0e-6*m23)) {
        *xc=(m12*x12-m23*x23-y12+y23)/(m12-m23);
        *yc=m12*(*xc-x12)+y12;
      }
      else { plerror("centre not found in routine tricen"); return FALSE; }
    }
    else { *xc=x23; *yc=m12*(x23-x12)+y12; }
  }
  else if (fabs(m23)>1e-6) {
    m23= -(x3-x2)/m23;
    *xc=x12; *yc=m23*(x12-x23)+y23;
  } else { plerror("centre not found in routine tricen"); return FALSE; }
  return TRUE;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *   P T B O X ( )      *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	ptbox(REAL *xc, *yc, short *i3, *j3, *k3)
 *
 *	Arguments:	(*xc, *yc)	= a given point which may need to be
 *				translated into the central periodic box
 *			*i3, *j3, *k3	= three periodic indices which we
 *				would like to be translated in tandem with
 *				the point (*xc, *yc)  (if need be).
 *
 *	Return value:	none
 *
 *	Action:		It is used by the Voronoi routines effectively to move
 *			a triangle into the central periodic box.
 *			The (*xc, *yc) are the given coordinates of the
 *			triangle's centre and *i3, *j3, *k3 are the periodic
 *			indices of the triangle vertices.  All of these given
 *			variables are changed so as to have (*xc, *yc) lie
 *			inside the central box.
 *
 *****************************************************************************/
void ptbox(xc, yc, i3, j3, k3)
REAL *xc, *yc;
short *i3, *j3, *k3;
{
  if (*xc< -boxwid/2.0) {
    *xc += boxwid;
    *i3=PERINCX(*i3); *j3=PERINCX(*j3); *k3=PERINCX(*k3);
  }
  else if (*xc>boxwid/2.0) {
    *xc -= boxwid;
    *i3=PERDECX(*i3); *j3=PERDECX(*j3); *k3=PERDECX(*k3);
  }
  if (*yc< -boxhgt/2.0) {
    *yc += boxhgt;
    *i3=PERINCY(*i3); *j3=PERINCY(*j3); *k3=PERINCY(*k3);
  }
  else if (*yc>boxhgt/2.0) {
    *yc -= boxhgt;
    *i3=PERDECY(*i3); *j3=PERDECY(*j3); *k3=PERDECY(*k3);
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *   I N B O X ( )      *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	boolean inbox(short i, i1)
 *
 *	Arguments:	i	= index of cell centre
 *				(= Delauney triangle apex)
 *			i1	= periodic index associated with i
 *
 *	Return value:	Tells whether Delauney triangle apex `i', after
 *			translating by the periodic vector encoded in `i1',
 *			lies **in** the central periodic **box**
 *
 *	Action:		(Ditto `Return value')
 *
 *****************************************************************************/
boolean inbox(i, i1)
short i, i1;
{
  REAL x, y;
  x=cx[i]+boxwid*PERX(i1); y=cy[i]+boxhgt*PERY(i1);
  if ((x<-boxwid/2.0) || (x>boxwid/2.0) || (y<-boxhgt/2.0) || (y>boxhgt/2.0)) return FALSE;
  else return TRUE;
}

/******************/
/* minor routines */
/******************/

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *     S Q F ( )        *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL sqf(REAL x)
 *
 *	Arguments:	x	= any real number
 *
 *	Return value:	x*x
 *
 *****************************************************************************/
REAL sqf(x)
REAL x; { return (x*x); }

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *    N O R M A L I Z E P ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	normalizep()
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		Shifts all the cell pressures and border pressures by
 *			a constant amount so that the average cell pressure
 *			works out at zero.
 *
 *			It is an arbitrary convention but it is useful because
 *			the (usually negative) values of Plateau border
 *			pressures tell us straightaway something about the
 *			curvature of Plateau border arcs.
 *			This routine is called whenever significant changes
 *			are made to the network in order to keep the pressures
 *			up to date.
 *
 *****************************************************************************/
void normalizep()
{
  /* This routine normalises the average cell pressure to zero */
  /* and sets 'bp' to have the appropriate negative value */
  short i, ii;
  REAL cptot=0.0, bptot=0.0, dbpav;
  for (ii=0; ii<nc; ii++) {
    i=clist[ii];
    cptot += cp[i];
  }
  cptot /= nc;
  /* now subtract this average pressure from everything! */
  for (ii=0; ii<nc; ii++) {
    i=clist[ii];
    cp[i] -= cptot;
  }
  for (ii=0; ii<nb; ii++) {
    i=blist[ii];
    bp[i] -= cptot;
    bptot += bp[i];
  }
  dbpav= -bpav;
  bpav=bptot/nb;
  dbpav += bpav;
  for (ii=0; ii<nb; ii++) {
    i=blist[ii];
    /* here relax the discrepancies from 'bpav' */
    bp[i] -= sbprelax*(bp[i]-bpav);
  }
  for (ii=0; ii<nbub; ii++) {
    i=bublist[ii];
    cp[i] += dbpav;
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *   N O R M A L I Z E B A ( )*     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	normalizeba()
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		Make sure that the total area of all the objects inside
 *			a periodic box is *equal* to the area of the periodic 
 *			box (discrepancies could arise due to cumulative
 *			numerical errors).  If an adjustment needs to be made
 *			then adjust the areas of the Plateau borders so that
 *			it all adds up properly.
 *
 *****************************************************************************/
void normalizeba()
{
  /* This subroutine ensures that the total area fits in the box by */
  /* fractionally adjusting the areas of the borders */
  short i, ii;
  REAL catot, batot, f;
  catot=0.0;
  for (ii=0; ii<nc; ii++) {
    i=clist[ii];
    catot += carea[i]+darea[i];
  }
  batot=0.0;
  for (ii=0; ii<nb; ii++) {
    i=blist[ii];
    batot += barea[i];
  }
  f=1.0+(boxwid*boxhgt-(catot+batot))/batot;
  for (ii=0; ii<nb; ii++) {
    i=blist[ii];
    barea[i] *= f;
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *P U T I N B O X ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	putinbox(short i)
 *
 *	Arguments:	i	= index of vertex in a wet foam
 *
 *	Return value:	none
 *
 *	Action:		Puts the vertex with index `i' back into the central
 *			periodic box (if it is not already in it).
 *			This ensures that the coordinates (vx[i], vy[i])
 *			always lie within the bounds of the central periodic
 *			box.  It is also necessary to update the array
 *			`vper[][]' at the same time (which tells you how to
 *			translate your neighbours so as to move them into the
 *			same vicinity as yourself).  The entries in `vper[][]'
 *			belonging to each of your three neighbours must also
 *			be updated (so that they still know how to move
 *			you into *their* vicinity).
 *
 *****************************************************************************/
void putinbox(i)
short i;
{
  REAL xlim, ylim, x, y;
  short j, k, k1;
  xlim=boxwid/2.0; ylim= boxhgt/2.0; 
  x=vx[i]; y=vy[i];
  if (x>xlim) {
    vx[i] -= boxwid;
    for (k=0; k<3; k++) {
      vper[i][k]=PERDECX(vper[i][k]);
      j=vnbr[i][k]; k1=knbr[k];
      /*
       * nbr moves in opposite direction...
       */
      vper[j][k1]=PERINCX(vper[j][k1]);
    }
  }
  else if (x<-xlim) {
    vx[i] += boxwid;
    for (k=0; k<3; k++) {
      vper[i][k]=PERINCX(vper[i][k]);
      j=vnbr[i][k]; k1=knbr[k];
      vper[j][k1]=PERDECX(vper[j][k1]);
    }
  }
  if (y>ylim) {
    vy[i] -= boxhgt;
    for (k=0; k<3; k++) {
      vper[i][k]=PERDECY(vper[i][k]);
      j=vnbr[i][k]; k1=knbr[k];
      vper[j][k1]=PERINCY(vper[j][k1]);
    }
  }
  else if (y<-ylim) {
    vy[i] += boxhgt;
    for (k=0; k<3; k++) {
      vper[i][k]=PERINCY(vper[i][k]);
      j=vnbr[i][k]; k1=knbr[k];
      vper[j][k1]=PERDECY(vper[j][k1]);
    }
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *   T R A N S ( )      *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	trans(REAL x, y, short per, REAL *x1, *y1)
 *
 *	Arguments:	(x, y)	= raw coordinate of vertex
 *			per	= periodic index (into which is coded a
 *				two-dimensional integer vector (per_x, per_y)
 *				which tells by how many periodic boxes the
 *				point (x, y) ought be translated)
 *			(*x1, *y1) = the resulting coordinates after being
 *				translated by the amount specified in the
 *				periodic index `per'.
 *
 *	Return value:	none
 *
 *	Action:		Translates the raw coordinate (x, y) through an
 *			integer number of periodic boxes up/down and left/right
 *			according to the integral vector encoded in the
 *			periodic index `per'.
 *
 *****************************************************************************/
void trans(x, y, per, x1, y1)
short per;
REAL x, y, *x1, *y1;
{
  *x1=x+boxwid*PERX(per);
  *y1=y+boxhgt*PERY(per);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  V N B R X Y ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	vnbrxy(short i, k, REAL *x1, *y1)
 *
 *	Arguments:	i	= index of a vertex in the wet foam
 *			k	= selects the k'th neighbour of vertex `i'
 *			(*x1, *y1) = the coordinates of the k'th nbr. of
 *				vertex `i', after having been translated into
 *				the immediate vicinity of vertex `i'.
 *
 *	Return value:	none
 *
 *	Action:		Return the coordinates of the k'th nbr. of vertex `i',
 *			after having moved the nbr. so as to lie in the
 *			immediate vicinity of vertex `i'.
 *
 *			NB NB NB -- This very simple routine is of pivotal
 *			importance in implementing periodic boundary conditions
 *			in the wet foam:  It is the most useful `black box'
 *			routine of all the periodic boundary routines.
 *			Whenever any geometrical calculations are done it is
 *			essential that all the points lie in the same vicinity
 *			(i.e. we shouldn't be using an image of the point
 *			which lies in the *wrong* periodic box).  This little
 *			subroutine guarantees that the neighbouring point
 *			(*x1, *y1) lies in the same vicinity as (vx[i], vy[i])
 *			and is therefore an essential first step in any
 *			geometrical calculation.
 *
 *****************************************************************************/
void vnbrxy(i, k, x1, y1)
short i, k;
REAL *x1, *y1;
{
  short j;
  j=vnbr[i][k];
  trans(vx[j], vy[j], vper[i][k], x1, y1);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *V O R T R A N S ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	vortrans(REAL x, y, xc, yc, *x1, *y1)
 *
 *	Arguments:	(x, y)	= raw coordinate which may need to be moved
 *			(xc, yc) = coordinate of point we want to be near to
 *			(*x1, *y1) = translated version of (x, y) after having
 *				been moved close to (xc, yc)
 *
 *	Return value:	none
 *
 *	Action:		A crude  method of translating the point (x, y) into
 *			the vicinity of (xc, yc) which has to be used in the
 *			Voronoi routine (in fact it is not guaranteed to work
 *			if the number of points is very small).  Essentially
 *			(x, y) is moved by an integer number of periodic boxes
 *			until it lies within half a periodic box of (xc, yc).
 *			The resulting point close to (xc, yc) is returned as
 *			the point (*x1, *y1).
 *
 *****************************************************************************/
void vortrans(x, y, xc, yc, x1, y1)
REAL x, y, xc, yc, *x1, *y1;
{
  REAL d, fsign();
  *x1=x;
  while (fabs(d= *x1-xc)>boxwid/2.0) *x1 -= fsign(d)*boxwid;
  *y1=y;
  while (fabs(d= *y1-yc)>boxhgt/2.0) *y1 -= fsign(d)*boxhgt;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *   F S I G N ( )      *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL fsign(REAL x)
 *
 *	Arguments:	x	= any real number
 *
 *	Return value:	x/|x|
 *
 *****************************************************************************/
REAL fsign(x)
REAL x;
{ return ((x>0) ? 1.0 : -1.0); }

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     V O R V N B R X Y ( )  *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	vorvnbrxy(short i, k, REAL *x1, *y1)
 *
 *	Arguments:	i	= index of vertex in a dry foam
 *			k	= selects the k'th nbr. of vertex `i'
 *			(*x1, *y1) = the coordinates of the k'th nbr. of
 *				vertex `i', after having been translated into
 *				the immediate vicinity of vertex `i'.
 *
 *	Return value:	none
 *
 *	Action:		Return the coordinates of the k'th nbr. of vertex `i',
 *			after having moved the nbr. so as to lie in the
 *			immediate vicinity of vertex `i'.
 *
 *			This is very similar to `vnbrxy()' except that it
 *			operates on a dry foam topology instead.
 *			(Note that this routine does *not* call `vortrans()'
 *			but calls `trans()' since topological information
 *			is available in array `vorvper[][]'.  Perhaps
 *			`vortrans()' would be better named ``crudetrans()'' ).
 *
 *****************************************************************************/
void vorvnbrxy(i, k, x1, y1)
short i, k;
REAL *x1, *y1;
{
  void trans();
  trans(vx[vorvnbr[i][k]], vy[vorvnbr[i][k]], vorvper[i][k], x1, y1);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     V O R K I N D E X ( )  *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	short vorkindex(short i, j, j1)
 *
 *	Arguments:	i	= index of first dry foam vertex
 *			j	= index of second dry foam vertex
 *			j1	= periodic index from `vorvper[i][*]' which
 *				would translate vertex `j' near to `i'
 *
 *	Return value:	an index `k' such that `i = vorvnbr[j][k]'
 *
 *	Action:		This routine essentially gives us the `inverse' of the
 *			array `vorvnbr[i][kk]'.  That is, suppose we know that
 *			`j' is the kk'th nbr. of `i' -- then this routine tells
 *			us that `i' is the k'th nbr. of `j', where `k' is the
 *			return value from the routine.
 *
 *			The apparently extraneous value `j1', giving the
 *			periodic index `vorvper[i][kk]', is necessary to cover
 *			a tricky eventuality.  It is theoretically possible
 *			that the *same* vertex `i' could be a neighbour of `j'
 *			more than once, since `j' could be connected to it in
 *			different periodic boxes.  Therefore to specify the
 *			neighbour relationship uniquely it is necessary to
 *			specify in which periodic box the nbr. lies via the
 *			periodic index `j1'.
 *
 *****************************************************************************/
short vorkindex(i, j, j1)
short i, j, j1;
{
  void plerror();
  short k1;
  j1=PERFN(-PERX(j1),-PERY(j1));
  for (k1=0; k1<3; k1++)
    if ((i==vorvnbr[j][k1]) && (j1==vorvper[j][k1])) return k1;
  plerror("no index found in routine vorkindex");
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *P E R I O D I C   I N D I C E S   * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutines:	perfn(short ix, iy), perx(short per), pery(short per),
 *			perincx(short per), perdecx(short per),
 *			perincy(short per), perdecy(short per).
 *
 *	Arguments:	perfn():  (ix, iy) = an integer vector representing
 *				translation through a certain number of
 8				periodic boxes.
 *			All others:    per = a periodic index
 *
 *	Return value:	none
 *
 *	Action:		Throughout the program it is frequently necessary to
 *			translate a vertex point (vx[], vx[]) through a
 *			certain number of periodic boxes.  This translation
 *			through a number of periodic boxes can be specified
 *			by a *periodic vector* (ix, iy) which is just a pair
 *			of integers each in the range [-8...+7].  However,
 *			these periodic vectors are not stored literally in
 *			arrays, instead they are combined and encoded into
 *			a *periodic index*.
 *
 *			A periodic index has information encoded in its bits
 *			in the following format:
 *			     _ _ _ _ _ _ _ _ _
 *			    |8|7|6|5|4|3|2|1|0|
 *			     - - - - - - - - -
 *			    | |       |       |
 *			     |     |       |
 *			     |     |       |
 *			     |     |       +--- x-vector component
 *			     |     +----------- y-vector component
 *			     +----------------- large-arc flag
 *			
 *			The ensuing subroutines simply make up a convenient
 *			interface to this compactly stored information and
 *			involve some simple exercises in binary arithmetic.
 *
 *			For example `perfn(ix, iy)' simply encodes the
 *			periodic vector (ix, iy) and returns a periodic index.
 *			The routines (perx(per), pery(per)) extract the
 *			periodic vector from the periodic index `per'.
 *			
 *			SETLARC() AND LARC():
 *			The routine `setlarc(per, la)' sets the large-arc flag
 *			bit of `per' to 0 or 1 according as `la' is true or
 *			false.  The complementary routine `larc(i, k)' returns
 *			the value (TRUE or FALSE) of the large-arc flag, for
 *			the arc which runs between vertex `i' and its k'th
 *			neighbour.
 *
 *****************************************************************************/
short perfn(ix,iy)
short ix,iy;
{
  return ( (iy & 0x0f)<<4 | (ix & 0x0f));
}

short perx(per)
short per;
{
  return ( (per&0x08) ? (per&0x0f)-16 : (per&0x0f) );
}

short pery(per)
short per;
{
  return ( (per&0x80) ? ((per&0xf0)>>4)-16 : (per&0xf0)>>4 );
}

short perincx(per)
short per;
{
  return ( (per & ~0x0f) | ((per & 0x0f)+1)%16 );
}

short perdecx(per)
short per;
{
  return ( (per & ~0x0f) | ((per & 0x0f)+15)%16 );
}

short perincy(per)
short per;
{
  return ( (per & ~0xf0) | (((per & 0xf0)+0x10) & 0xf0) );
}

short perdecy(per)
short per;
{
  return ( (per & ~0xf0) | (((per & 0xf0)+0xf0) & 0xf0) );
}

boolean larc(i, k)
short i, k;
{
  return (boolean) (vper[i][k] & LARC)>>8;
}
  
void setlarc(i1, la)
short *i1;
boolean la;
{ *i1=(*i1 & PERMASK) | ((la) ? LARC : 0x00); }

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    L I S T S   O F   V A L I D   O B J E C T S  * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutines:	short bgetindex(), short vgetindex(),
 * 			vforgetindex(short i), cforgetindex(short i),
 *			bforgetindex(short i).
 *
 *	Arguments:	i	= index of object to be forgotten.
 *
 *	Return value:	none
 *
 *	Action:		Since the topology of the network changes throughout
 *			its life, the number of objects of different types
 *			(cells, vertices etc.) will also change.  The lists
 *			of valid indices for various objects are maintained
 *			in the following arrays:
 *				vlist[0..nv] = list of valid vertex indices
 *				clist[0..nc] = list of valid cell indices
 *				blist[0..nb] = list of valid border indices
 *			The indices *stored* within these arrays lie within
 *			the ranges [0..onv], [0..onc] and [0..onb]
 *			respectively, corresponding to the original (maximum)
 *			sizes of these arrays.
 *
 *			The following few subroutines manage the lists of
 *			valid indices.  The `*forgetindex(i)' routines delete
 *			index `i' from list `*list[]' and the `*getindex()'
 *			routines return an unused index from list `*list[]'
 *			for use as a new object.
 *
 *****************************************************************************/
short vgetindex()
{
  void plerror();
  short ii, i, j, kk;
  j= -1;
  if (nv<onv) {
    for (ii=0; ii<nv+1; ii++) {
      i=vlist[ii];
      if (i>j+1 || ii==nv) {
        for (kk=nv; kk>ii; kk--) vlist[kk]=vlist[kk-1];
        vlist[ii]=j+1; nv++;
        return(j+1);
      }
      j=i;
    }
  }
  plerror("no index found in routine vgetindex");
}

short bgetindex()
{
  void plerror();
  short i, j, ii, kk;
  j= -1;
  if (nb<onb) {
    for (ii=0; ii<nb+1; ii++) {
      i=blist[ii];
      if (i>j+1 || ii==nb) {
        for (kk=nb; kk>ii; kk--) blist[kk]=blist[kk-1];
        blist[ii]=j+1; nb++;
        return(j+1);
      }
      j=i;
    }
  }
  plerror("no index found in routine bgetindex");
}

void vforgetindex(i)
short i;
{
  void plerror();
  short ii, kk;
  for (ii=0; ii<nv; ii++)
    if (i==vlist[ii]) break;
  if (ii==nv) 
    plerror("index not found in routine vforgetindex()");
  else {
    for (kk=ii; kk<nv-1; kk++) vlist[kk]=vlist[kk+1];
    nv--;
    vlist[nv]=onv+1;
  }
}

void cforgetindex(i)
short i;
{
  void plerror();
  short ii, kk;
  for (ii=0; ii<nc; ii++)
    if (i==clist[ii]) break;
  if (ii==nc) 
    plerror("index not found in routine cforgetindex()");
  else {
    for (kk=ii; kk<nc-1; kk++) clist[kk]=clist[kk+1];
    nc--;
    clist[nc]=onc+1;
  }
}

void bforgetindex(i)
short i;
{
  void plerror();
  short ii, kk;
  for (ii=0; ii<nb; ii++)
    if (i==blist[ii]) break;
  if (ii==nb) 
    plerror("index not found in routine bforgetindex()");
  else {
    for (kk=ii; kk<nb-1; kk++) blist[kk]=blist[kk+1];
    nb--;
    blist[nb]=onb+1;
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *B U B I N B O X ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	bubinbox(REAL *xc, *yc)
 *
 *	Arguments:	(*xc, *yc)	= centre of isolated bubble
 *
 *	Return value:	none
 *
 *	Action:		Translates the centre point of an isolated bubble into
 *			the central periodic box, if it is not already in
 *			there.
 *
 *****************************************************************************/
void bubinbox(xc, yc)
REAL *xc, *yc;
{
  REAL halfwid, halfhgt;
  REAL fsign();
  halfwid=boxwid/2.0; halfhgt=boxhgt/2.0;
  while (fabs(*xc)>halfwid) *xc -= fsign(*xc)*boxwid;
  while (fabs(*yc)>halfhgt) *yc -= fsign(*yc)*boxhgt;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     P E R C O N C A T ( )  *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	short perconcat(short i1, j1)
 *
 *	Arguments:	i1, j1	= two periodic indices
 *
 *	Return value:	a `concatenated' periodic index
 *
 *	Action:		This routine adds the periodic vectors encoded in `i1'
 *			and `j1', then returns the sum encoded in a new
 *			periodic vector. (Note: any information stored in the
 *			`large-arc' flag bits of `i1' and `j1' is not carried
 *			over to the value returned).
 *
 *****************************************************************************/
short perconcat(i1, j1)
short i1, j1;
{
  short ix, iy, jx, jy;
  ix=PERX(i1); iy=PERY(i1);
  jx=PERX(j1); jy=PERY(j1);
  return PERFN(ix+jx, iy+jy);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *   P I M O D ( )      *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL pimod(REAL a)
 *
 *	Arguments:	a	= given angle in radians
 *
 *	Return value:	angle in range [-PI,PI]
 *
 *	Action:		Takes the given angle `a' and maps it into the
 *			range [-PI,PI].
 *
 *****************************************************************************/
REAL pimod(a)
REAL a;
{
  while (fabs(a)>PI) a -= 2.0*PI*( (a>0.0) ? 1.0 : -1.0);
  return a;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     T W O P I M O D ( )    *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL twopimod(REAL a)
 *
 *	Arguments:	a	= given angle in radians
 *
 *	Return value:	angle in range [0,2*PI]
 *
 *	Action:		Takes the given angle `a' and maps it into the
 *			range [0,2*PI].
 *
 *****************************************************************************/
REAL twopimod(a)
REAL a;
{
  boolean lt=FALSE, gt=FALSE;
  REAL pi2;
  pi2=2.0*PI;
  do {
    a += pi2*(( (lt=(a<0.0)) ? 1.0 : 0.0) - ((gt=(a>pi2)) ? 1.0 : 0.0));
  } while (lt || gt);
  return a;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  L I N L E N ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL linlen(REAL x1, y1, x2, y2)
 *
 *	Arguments:	(x1, y1), (x2, y2) = endpoints of a line segment
 *
 *	Return value:	length of the line segment
 *
 *****************************************************************************/
REAL linlen(x1, y1, x2, y2)
REAL x1, y1, x2, y2;
{
  return (sqrt(((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))) );
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  L I N L E N 2 ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL linlen2(REAL x1, y1, x2, y2)
 *
 *	Arguments:	(x1, y1), (x2, y2) = endpoints of a line segment
 *
 *	Return value:	the square of the length of the line segment
 *
 *****************************************************************************/
REAL linlen2(x1, y1, x2, y2)
REAL x1, y1, x2, y2;
{
  return((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     C A R C A N G L E ( )  *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL carcangle(REAL x1, y1, x2, y2, p1, p2)
 *
 *	Arguments:	(x1, y1)	= first endpoint of cell-cell arc
 *			(x2, y2)	= second endpoint of cell-cell arc
 *			p1,  p2		= the pressures on either side of the
 *				cell-cell arc
 *
 *	Return value:	The angle subtended by the arc joining point (x1, y1)
 *			to (x2, y2), where the curvature of the arc is
 *			determined by the difference in cell pressures p2 - p1.
 *			Sign convention for returned angle `theta':
 *				p1 > p2  ==>  theta > 0
 *				p1 < p2  ==>  theta < 0
 *
 *****************************************************************************/
REAL carcangle(x1, y1, x2, y2, p1, p2)
REAL x1, y1, x2, y2, p1, p2;
{
  void plerror();
  REAL dp, x, linlen();
  if (fabs(dp=p1-p2)<MPRECISION)
    return 0.0;
  else
    if (fabs(x=(0.5/CRADIUS(p1,p2))*linlen(x1, y1, x2, y2))>1.0)
      plerror("arcbreakage in routine carcangle");
    else return(2.0*asin(x));
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     C A R C L E N ( )      *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL carclen(REAL x1, y1, x2, y2, p1, p2)
 *
 *	Arguments:	(x1, y1)	= first endpoint of cell-cell arc
 *			(x2, y2)	= second endpoint of cell-cell arc
 *			p1,  p2		= the pressures on either side of the
 *				cell-cell arc
 *
 *	Return value:	The arc length  of the cell-cell arc joining (x1, y1)
 *			to (x2, y2), where the curvature of the arc is
 *			determined by the difference in cell pressures p2 - p1.
 *
 *****************************************************************************/
REAL carclen(x1, y1, x2, y2, p1, p2)
REAL x1, y1, x2, y2, p1, p2;
{
  REAL dp, linlen(), carcangle();
  if (fabs(dp=p1-p2)<MPRECISION)
    return linlen(x1, y1, x2, y2);
  else
    return(carcangle(x1, y1, x2, y2, p1, p2)*CRADIUS(p1,p2) );
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     B A R C A N G L E ( )  *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL barcangle(REAL x1, y1, x2, y2, p1, p2,
 *				boolean la, *arcbrk)
 *
 *	Arguments:	(x1, y1)	= first endpoint of Plateau border arc
 *			(x2, y2)	= second endpoint of Plateau border arc
 *			p1		= pressure of adjacent cell
 *			p2		= pressure of adjacent Plateau border
 *			la	= large-arc flag
 *			*arcbrk = flag to tell if `arc breakage' occurred
 *
 *	Return value:	The angle subtended by the arc joining point (x1, y1)
 *			to (x2, y2), where the curvature of the arc is
 *			determined by the difference between cell pressure `p1'
 *			and Plateau border pressure `p2'.
 *
 *			LARGE AND SMALL ARCS:
 *			There is an ambiguity in the correct value to return
 *			for this angle.  Given the radius of the border arc,
 *			determined by the macro BRADIUS(p1,p2), there are *two*
 *			possible arcs which could be fitted between the end
 *			points.  Either a small arc with `theta < PI' or a
 *			large arc `2*PI - theta > PI'.  When the large arc flag
 *			`la == TRUE' then the larger angle is selected.
 *
 *			`ARC BREAKAGE' ERROR:
 *			An important error condition occurs if 2*BRADIUS(p1,p2)
 *			is smaller than the distance between the endpoints
 *			(x1, y1) and (x2, y2).  In that case it is impossible
 *			to fit an arc between the endpoints and the error
 *			condition `arc breakage' is signalled by setting
 *			`*arcbrk = TRUE'.  THIS ERROR MUST ALWAYS BE CHECKED.
 *			If you ignore it then the network is guaranteed to
 *			crash (very rapidly).
 *
 *****************************************************************************/
REAL barcangle(x1, y1, x2, y2, p1, p2, la, arcbrk)
REAL x1, y1, x2, y2, p1, p2;
boolean la, *arcbrk;
{
  REAL dp, x, linlen(), fsign();
  void plerror();
  *arcbrk=FALSE;
  if (fabs(dp=p1-p2)<MPRECISION)
    return((la) ? 2.0*PI : 0.0);
  else
    if (fabs(x=(0.5/BRADIUS(p1,p2))*linlen(x1,y1,x2,y2))>1.0) {
      *arcbrk=TRUE;
#ifndef BACKGROUND
      plerror("arcbreakage in routine barcangle");
#endif
    }
    else {
      x=2.0*asin(x);
      return((la) ? fsign(x)*(2.0*PI-fabs(x)) : x);
    }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     * B A R C L E N ( )    *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL barclen(REAL x1, y1, x2, y2, p1, p2,
 *				boolean la, *arcbrk)
 *
 *	Arguments:	(x1, y1)	= first endpoint of Plateau border arc
 *			(x2, y2)	= second endpoint of Plateau border arc
 *			p1		= pressure of adjacent cell
 *			p2		= pressure of adjacent Plateau border
 *			la	= large-arc flag
 *			*arcbrk = flag to tell if `arc breakage' occurred
 *
 *	Return value:	Returns the arc length of the arc joining (x1, y1) to
 *			(x2, y2) whose curvature is determined by the
 *			difference of pressures `p1' and `p2'.  See routine
 *			`barcangle()' above for comments on the significance
 *			of `la',  `*arcbrk' and an important error condition.
 *
 *****************************************************************************/
REAL barclen(x1, y1, x2, y2, p1, p2, la, arcbrk)
REAL x1, y1, x2, y2, p1, p2;
boolean la, *arcbrk;
{
  REAL dp, linlen(), barcangle();
  *arcbrk=FALSE;
  if (fabs(dp=p1-p2)<MPRECISION)
    return linlen(x1,y1,x2,y2);
  else
    return(barcangle(x1,y1,x2,y2,p1,p2,la,arcbrk)*BRADIUS(p1,p2));
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *     G E E ( )        *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL gee(REAL alpha)
 *
 *	Arguments:	alpha	= angle subtending an arc
 *
 *	Return value:	Gives the ratio of the arc length to the length of the
 *			chord of the arc.
 *
 *****************************************************************************/
REAL gee(alpha)
REAL alpha;
{
  return(0.5*alpha/sin(0.5*alpha));
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     A R C C E N T R E ( )  *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	boolean arccentre(x1, y1, x2, y2, alpha, *xc, *yc)
 *
 *	Arguments:	(x1, y1)	= first endpoint of arc
 *			(x2, y2)	= second endpoint of arc
 *			alpha		= angle subtending the arc
 *			(*xc, *yc)	= centre point of the arc
 *
 *	Return value:	Returns FALSE if the angle `alpha' is small and, in
 *			that case, does not put the centre point coordinates
 *			in (*xc, *yc) but instead sets them to the *midpoint*
 *			of the line (x1, y1), (x2, y2).
 *
 *	Action:		Either puts the coordinate of the arc centre into
 *			the coordinates (*xc, *yc) and returns TRUE (when
 *			`alpha' is reasonably large)
 *					or
 *			puts the coodinates of the midpoint of the line (x1,y1)
 *			to (x2, y2) into the variabls (*xc, *yc) and returns
 *			the value FALSE (when `alpha' is fairly small).
 *
 *****************************************************************************/
boolean arccentre(x1, y1, x2, y2, alpha, xc, yc)
REAL x1, y1, x2, y2, alpha, *xc, *yc;
{
  REAL dx, dy, t;
  if (fabs(alpha)<0.05) {
    *xc=(x1+x2)/2.0; *yc=(y1+y2)/2.0;
    return FALSE;
  }
  else {
    dx=0.5*(x2-x1); dy=0.5*(y2-y1);
    t=tan(0.5*alpha);
    *xc=x1+dx+dy/t; *yc=y1+dy-dx/t;
    return TRUE;
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     * T R I A R E A ( )    *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL triarea(REAL x1, y1, x2, y2, x3, y3)
 *
 *	Arguments:	(x1, y1)	= first apex
 *			(x2, y2)	= second apex
 *			(x3, y3)	= third apex
 *
 *	Return value:	Area of triangle with given apices.
 *
 *			Sign convention: Positive when the three apices are
 *			in anti-clockwise order, else negative.
 *
 *****************************************************************************/
REAL triarea(x1, y1, x2, y2, x3, y3)
REAL x1, y1, x2, y2, x3, y3;
{
  x2 -= x1; y2 -= y1; x3 -= x1; y3 -= y1;
  return(0.5*(x2*y3-x3*y2));
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     Q U A D A R E A ( )    *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL quadarea(REAL x1, y1, x2, y2, x3, y3, x4, y4)
 *
 *	Arguments:	(x1, y1)	= first apex
 *			(x2, y2)	= second apex
 *			(x3, y3)	= third apex
 *			(x4, y4)	= fourth apex
 *
 *	Return value:	Area of quadrangle with given vertices.
 *
 *			Sign convention: Positive when the four vertices are
 *			in anti-clockwise order.
 *
 *****************************************************************************/
REAL quadarea(x1, y1, x2, y2, x3, y3, x4, y4)
REAL x1, y1, x2, y2, x3, y3, x4, y4;
{
  REAL triarea();
  return(triarea(x1,y1,x2,y2,x3,y3)+triarea(x1,y1,x3,y3,x4,y4));
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *C A R C A R E A ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL carcarea(REAL alpha, p1, p2)
 *
 *	Arguments:	alpha	= the angle subtending the cell-cell arc
 *			p1, p2	= the pressures on either side of the arc
 *
 *	Return value:	Returns the area between an arc and its chord.  The
 *			area returned has the same sign as `alpha'
 *			(irrespective of the order in which `p1' and `p2' are
 *			passed).
 *
 *****************************************************************************/
REAL carcarea(alpha, p1, p2)
REAL alpha, p1, p2;
{
  REAL dp, r;
  if (fabs(dp=p1-p2)<MPRECISION)
    return 0.0;
  else {
    r=CRADIUS(p1,p2);
    return(r*r*(0.5*alpha-sin(0.5*alpha)*cos(0.5*alpha)));
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *B A R C A R E A ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL barcarea(REAL alpha, p1, p2)
 *
 *	Arguments:	alpha	= the angle subtending the border arc
 *			p1, p2	= the pressures on either side of the arc
 *
 *	Return value:	Returns the area between an arc and its chord.  The
 *			area returned has the same sign as `alpha'
 *			(irrespective of the order in which `p1' and `p2' are
 *			passed).
 *
 *****************************************************************************/
REAL barcarea(alpha, p1, p2)
REAL alpha, p1, p2;
{
  REAL dp, r;
  if (fabs(dp=p1-p2)<MPRECISION)
    return 0.0;
  else {
    r=BRADIUS(p1,p2);
    return(r*r*(0.5*alpha-sin(0.5*alpha)*cos(0.5*alpha)));
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *L I N A N G L E ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL linangle(REAL x1, y1, x2, y2)
 *
 *	Arguments:	(x1, y1),  (x2, y2) = endpoints of line
 *
 *	Return value:	Returns the angle between the given line and the
 *			positive x-axis (in the range [0..2*PI] ).
 *
 *****************************************************************************/
REAL linangle(x1, y1, x2, y2)
REAL x1, y1, x2, y2;
{
  REAL twopimod();
  if (fabs(x2-x1)<MPRECISION)
    return twopimod(0.5*PI*fsign(y2-y1));
  else
    return twopimod(atan((y2-y1)/(x2-x1))+((x2<x1) ? PI : 0.0));
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  C A N G L E ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL cangle(REAL x1, y1, x2, y2, p1, p2)
 *
 *	Arguments:	
 *				    p2
 *				  -------
 *				/        \
 *			  (x1, y1)      (x2, y2)
 *				    p1
 *
 *			(x1, y1), (x2, y2) = endpoints of the cell-cell arc
 *			p1, p2	= pressures of the adjacent cells as
 *				illustrated roughly in the diagram
 *
 *	Return value:	Gives the angle between the positive x-axis and the
 *			tangent to the arc at (x1, y1).  
 *
 *****************************************************************************/
REAL cangle(x1, y1, x2, y2, p1, p2)
REAL x1, y1, x2, y2, p1, p2;
{
  REAL carcangle(), twopimod(), x;
  x=twopimod(linangle(x1,y1,x2,y2)+0.5*carcangle(x1,y1,x2,y2,p1,p2));
  return x;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  B A N G L E ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL bangle(REAL x1, y1, x2, y2, p1, p2,
 *				boolean la, *arcbrk)
 *
 *	Arguments:	
 *				    p2
 *				  -------
 *				/        \
 *			  (x1, y1)      (x2, y2)
 *				    p1
 *
 *			(x1, y1), (x2, y2) = endpoints of the border arc
 *			p1, p2	= pressures of adjacent cell/Plateau border
 *				as illustrated roughly in the diagram
 *			la	= large-arc flag
 *			*arcbrk = flag to indicate `arc breakage' condition
 *
 *	Return value:	Gives the angle between the positive x-axis and the
 *			tangent to the arc at (x1, y1).  
 *
 *			See routine `barcangle()' above for comments on the
 *			significance of arguments `la' and `*arcbrk'.
 *
 *****************************************************************************/
REAL bangle(x1, y1, x2, y2, p1, p2, la, arcbrk)
REAL x1, y1, x2, y2, p1, p2;
boolean la, *arcbrk;
{
  REAL barcangle(), twopimod(), x;
  *arcbrk=FALSE;
  x=barcangle(x1,y1,x2,y2,p1,p2,la,arcbrk);
  if (! *arcbrk) {
    x=twopimod(linangle(x1,y1,x2,y2)+0.5*x);
    return x;
  }
  else {
    *arcbrk=TRUE;
    return 0.0;
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *B C A N G L E 1 ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL bcangle1(REAL x1, y1, x2, y2, x3, y3, p1, p2, pb,
				 boolean la, arcbrk)
 *
 *	Arguments:	(x1, y1)	= first point of border arc
 *			(x2, y2)	= second point of border arc
 *					(= first point of cell-cell arc)
 *			(x3, y3)	= second point of cell-cell arc
 *			p1	= pressure of the cell which lies on the
 *				*same* side of the arc (x2, y2), (x3, y3)
 *				as the point (x1, y1) does
 *			p2	= pressure of the cell on the *other* side
 *				of the arc (x2, y2), (x3, y3)
 *			pb	= pressure of the Plateau border which lies
 *				on the opposite side of the border arc
 *				(x1, y1), (x2, y2) to the cell with
 *				pressure `p1'.
 *			la	= large-arc flag for arc (x1, y1), (x2, y2)
 *			*arcbrk = `arc breakage' flag for the border arc
 *
 *
 *	Return value:	Returns the difference between the border and cell-cell
 *			arc tangent angles, at the point (x2, y2).  When the
 *			slopes of these two arcs match smoothly at this point
 *			then this routine would return exactly PI.  This
 *			routine covers one orientation of the three points --
 *			the other orientation is covered by `bcangle2()'.
 *			(e.g. if (x1, y1), (x2, y2), (x3, y3) are taken to be
 *			the apices of a triangle then this routine assumes
 *			that they lie in clockwise order).
 *
 *			See routine `barcangle()' above for comments on the
 *			significance of arguments `la' and `*arcbrk'.
 *
 *****************************************************************************/
REAL bcangle1(x1, y1, x2, y2, x3, y3, p1, p2, pb, la, arcbrk)
REAL x1, y1, x2, y2, x3, y3, p1, p2, pb;
boolean la, *arcbrk;
{
  REAL bangle(), cangle(), twopimod(), x;
  x=cangle(x2,y2,x3,y3,p1,p2)-bangle(x2,y2,x1,y1,pb,p1,la,arcbrk);
  return( (*arcbrk) ? 0.0 : twopimod(x) );
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *B C A N G L E 2 ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL bcangle2(REAL x1, y1, x2, y2, x3, y3, p1, p2, pb,
				 boolean la, arcbrk)
 *
 *	Arguments:	(x1, y1)	= first point of border arc
 *			(x2, y2)	= second point of border arc
 *					(= first point of cell-cell arc)
 *			(x3, y3)	= second point of cell-cell arc
 *			p1	= pressure of the cell which lies on the
 *				*same* side of the arc (x2, y2), (x3, y3)
 *				as the point (x1, y1) does
 *			p2	= pressure of the cell on the *other* side
 *				of the arc (x2, y2), (x3, y3)
 *			pb	= pressure of the Plateau border which lies
 *				on the opposite side of the border arc
 *				(x1, y1), (x2, y2) to the cell with
 *				pressure `p1'.
 *			la	= large-arc flag for arc (x1, y1), (x2, y2)
 *			*arcbrk = `arc breakage' flag for the border arc
 *
 *
 *	Return value:	Returns the difference between the border and cell-cell
 *			arc tangent angles, at the point (x2, y2).  When the
 *			slopes of these two arcs match smoothly at this point
 *			then this routine would return exactly PI.  This
 *			routine covers one orientation of the three points --
 *			the other orientation is covered by `bcangle1()'.
 *			(e.g. if (x1, y1), (x2, y2), (x3, y3) are taken to be
 *			the apices of a triangle then this routine assumes
 *			that they lie in anti-clockwise order).
 *
 *			See routine `barcangle()' above for comments on the
 *			significance of arguments `la' and `*arcbrk'.
 *
 *****************************************************************************/
REAL bcangle2(x1, y1, x2, y2, x3, y3, p1, p2, pb, la, arcbrk)
REAL x1, y1, x2, y2, x3, y3, p1, p2, pb;
boolean la, *arcbrk;
{
  REAL bangle(), cangle(), twopimod(), x;
  x=cangle(x2,y2,x3,y3,p2,p1)-bangle(x2,y2,x1,y1,p1,pb,la,arcbrk);
  return( (*arcbrk) ? 0.0 : twopimod(x) );
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *C E L L A R E A ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL cellarea(short i, k, boolean *arcbrk)
 *
 *	Arguments:	i	= index of vertex lying on perimeter of cell
 *			k	= select k'th cell adjacent to vertex `i'
 *			*arcbrk = `arc breakage' flag
 *
 *	Return value:	Area of the cell specified by `i' and `k'.
 *			The regions adjacent to vertex `i' are given in a
 *			standard order by the neighbour index `k = 0, 1 or 2'.  
 *			The value `k=0' always selects the adjacent Plateau
 *			border (an illegal choice for this routine) and `k=1'
 *			or `k=2' select the two adjacent cells in
 *			anti-clockwise order about vertex `i'.
 *
 *			See routine `barcangle()' above for comments on the
 *			significance of argument `*arcbrk'.
 *
 *****************************************************************************/
REAL cellarea(i,k,arcbrk)
short i, k;
boolean *arcbrk;
{
  REAL ca, cl;
  void careaperim();
  careaperim(i,k,TRUE,FALSE,&ca,&cl,arcbrk);
  return ca;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *C E L L P E R I M ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL cellperim(short i, k, boolean *arcbrk)
 *
 *	Arguments:	i	= index of vertex lying on perimeter of cell
 *			k	= select k'th cell adjacent to vertex `i'
 *			*arcbrk = `arc breakage' flag
 *
 *	Return value:	Perimeter length of the cell specified by `i' and `k'.
 *			The regions adjacent to vertex `i' are given in a
 *			standard order by the neighbour index `k = 0, 1 or 2'.  
 *			The value `k=0' always selects the adjacent Plateau
 *			border (an illegal choice for this routine) and `k=1'
 *			or `k=2' select the two adjacent cells in
 *			anti-clockwise order about vertex `i'.
 *
 *			See routine `barcangle()' above for comments on the
 *			significance of argument `*arcbrk'.
 *
 *****************************************************************************/
REAL cellperim(i,k,arcbrk)
short i, k;
boolean *arcbrk;
{
  REAL ca, cl;
  void careaperim();
  careaperim(i,k,FALSE,TRUE,&ca,&cl,arcbrk);
  return cl;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *    C E L L A R E A P E R I M ( ) * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	careaperim(short i, k, boolean caflg, clflg,
 *				REAL *ca, *cl, boolean *arcbrk)
 *
 *	Arguments:	i	= index of vertex lying on perimeter of cell
 *			k	= select k'th cell adjacent to vertex `i'
 *			caflg	= flag to calculate cell area
 *			clflg	= flag to calculate cell perimeter length
 *			*ca	= returns cell area (if `caflg' was TRUE)
 *			*cl	= returns cell perimeter (if `clflg' was TRUE)
 *			*arcbrk = `arc breakage' flag
 *
 *	Return value:	Calculates, optionally, both the cell area and
 *			perimeter length for the cell specified by `i' and `k'.
 *			The regions adjacent to vertex `i' are given in a
 *			standard order by the neighbour index `k = 0, 1 or 2'.  
 *			The value `k=0' always selects the adjacent Plateau
 *			border (an illegal choice for this routine) and `k=1'
 *			or `k=2' select the two adjacent cells in
 *			anti-clockwise order about vertex `i'.
 *
 *			See routine `barcangle()' above for comments on the
 *			significance of argument `*arcbrk'.
 *
 *****************************************************************************/
void careaperim(i, k, caflg, clflg, ca, cl, arcbrk)
REAL *ca, *cl;
short i, k;
boolean caflg, clflg, *arcbrk;
{
  short i1, ii, ii1, j, j1, c, perconcat();
  REAL tx, ty, p1, p2, x1, y1, x2, y2, x3, y3, alpha, pb,
    barcangle(), carcangle(), barclen(), barcarea(), carcarea(), triarea(),
    fsign();
  boolean la, larc();
  void trans();
  *arcbrk=FALSE;
  if (k==2) { i=vnbr[i][0]; k=1; }
  c=cadj[i][k];
  ii=i; p1=cp[c]; *ca=0.0; *cl=0.0;
  i1=0; j=vnbr[ii][2]; j1=vper[ii][2];
  tx=(REAL) PERX(j1); ty=(REAL) PERY(j1);
  x1=vx[ii]; y1=vy[ii]; x2=x1; y2=y1;
  do {
    x3=vx[j]+tx*boxwid; y3=vy[j]+ty*boxhgt;
    la=larc(i,2); pb=bp[cadj[i][0]];
    alpha=barcangle(x3,y3,x2,y2,p1,pb,la,arcbrk);
    if (*arcbrk) break;
    if (caflg) *ca += barcarea(alpha,p1,pb);
    if (clflg) *cl += barclen(x3,y3,x2,y2,p1,pb,la,arcbrk);
    if (*arcbrk) break;
    if ((i!=ii) && caflg) *ca += triarea(x1,y1,x2,y2,x3,y3);
    i=j; i1=j1;
    tx += PERX(vper[j][0]); ty += PERY(vper[j][0]);
    j=vnbr[j][0];
    x2=x3; y2=y3;
    x3=vx[j]+tx*boxwid; y3=vy[j]+ty*boxhgt;
    if ((j!=ii) && caflg) *ca += triarea(x1,y1,x2,y2,x3,y3);
    p2=cp[cadj[i][1]];
    alpha=carcangle(x3,y3,x2,y2,p1,p2);
    if (caflg) *ca += carcarea(alpha,p1,p2);
    if (clflg) *cl += carclen(x3,y3,x2,y2,p1,p2);
    if (j==ii) break;
    i=j; i1=j1;
    tx += PERX(vper[j][2]); ty += PERY(vper[j][2]);
    j=vnbr[j][2];
    x2=x3; y2=y3;
  } while (TRUE);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *B O R D A R E A ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL bordarea(short i, boolean *arcbrk)
 *
 *	Arguments:	i	= index of vertex lying on Plateau border
 *			*arcbrk = `arc breakage' flag
 *
 *	Return value:	Returns the area of the Plateau border that vertex `i'
 *			is adjacent to (note that unlike in the case of cells,
 *			this is sufficient to specify the border uniquely).
 *
 *			See routine `barcangle()' above for comments on the
 *			significance of argument `*arcbrk'.
 *
 *****************************************************************************/
REAL bordarea(i,arcbrk)
short i;
boolean *arcbrk;
{
  short i1, ii, ii1, j, j1, b, perconcat();
  REAL ba, p1, p2, x1, y1, x2, y2, x3, y3, alpha, barcangle(),
    barclen(), barcarea(), triarea(), fsign();
  boolean la, larc();
  void trans();
  *arcbrk=FALSE;
  b=cadj[i][0];
  ii=i; p1=bp[b];
  i1=0; j=vnbr[ii][1]; j1=vper[ii][1];
  x1=vx[ii]; y1=vy[ii]; x2=x1; y2=y1;
  ba=0.0;
  do {
    trans(vx[j],vy[j],j1,&x3,&y3);
    la=larc(i,1); p2=cp[cadj[i][2]];
    alpha=barcangle(x3,y3,x2,y2,p1,p2,la,arcbrk);
    if (*arcbrk) break;
    ba += barcarea(alpha,p2,p1);
    if (j==ii) break;
    if (i!=ii) ba += triarea(x1,y1,x2,y2,x3,y3);
    i=j; i1=j1;
    j1=perconcat(j1,vper[j][1]);
    j=vnbr[j][1];
    x2=x3; y2=y3;
  } while (TRUE);
  return (ba);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *    P E R C U L A T E S ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	boolean perculates(short i)
 *
 *	Arguments:	i	= index of a vertex lying on a Plateau border
 *
 *	Return value:	TRUE if percolation were to occur once the cell-cell
 *			arc adjacent to vertex `i' were removed.
 *
 *	Action:		This routine detects a Plateau border on the brink
 *			of percolating across the network, using a trick.  It
 *			relies on the fact that if you jump from vertex to
 *			vertex all along the perimeter of a Plateau border, you
 *			usually end up where you started.  However, if you do
 *			this with a *percolating* Plateau border then, when you
 *			get back to the first vertex, you will find that you
 *			have travelled one periodic box away from where you
 *			started !  Therefore by keeping track of the periodic
 *			vector accumulated as you traverse the perimeter,
 *			you can detect if the Plateau border is percolating.
 *
 *****************************************************************************/
boolean perculates(i)
short i;
{
  /* This routine answers the question: if the edge adjacent to 'i' were */
  /* removed, would the resulting PB percolate? */
  short j, j1;
  short perconcat();
  j=i; j1=0;
  do {
    j1=perconcat(j1,vper[j][1]); j=vnbr[j][1];
  } while (vnbr[j][0]!=i);
  return(j1 != vper[i][0]);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *    D C A D J A R E A S ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	dcadjareas(short i, REAL dbp)
 *
 *	Arguments:	i	= index of vector lying on Plateau border
			dbp	= increment to Plateau border pressure
 *
 *	Return value:	none
 *
 *	Action:		Adjusts the recorded geometrical sizes of all the
 *			cells adjacent to the Plateau border (in `carea[c]'),
 *			assuming that this individual Plateau border has just
 *			had its pressure raised by the amount `dbp'.
 *			Note that the *target* areas of the cells are not
 *			changed (that is, `carea[c]+darea[c] = unchanged')
 *			and this routine does not actually carry out
 *			any changes to the border pressures `bp[]'.
 *
 *****************************************************************************/
void dcadjareas(i,dbp)
short i;
REAL dbp;
{
  short i1, ii, ii1, j, j1, b, c, perconcat();
  REAL dca, p1, p2, x1, y1, x2, y2, x3, y3, alpha1, alpha2,
    barcangle(), barclen(), barcarea(), triarea(), fsign();
  boolean la, larc(), arcbrk;
  void trans(), cellpop();
  b=cadj[i][0];
  ii=i; p2=bp[b];
  i1=0; j=vnbr[ii][1]; j1=vper[ii][1];
  x1=vx[ii]; y1=vy[ii]; x2=x1; y2=y1;
  do {
    trans(vx[j],vy[j],j1,&x3,&y3);
    la=larc(i,1); p1=cp[c=cadj[i][2]];
    alpha1=barcangle(x3,y3,x2,y2,p1,p2,la,&arcbrk);
    if (arcbrk) {
      cellpop(i,2,pressuredelta/2.0);
      alpha1=barcangle(x3,y3,x2,y2,p1,p2,la,&arcbrk);
    }
    alpha2=barcangle(x3,y3,x2,y2,p1,p2+dbp,la,&arcbrk);
    if (arcbrk) {
      cellpop(i,2,pressuredelta/2.0);
      alpha2=barcangle(x3,y3,x2,y2,p1,p2+dbp,la,&arcbrk);
    }
    dca=barcarea(alpha2,p2,p1)-barcarea(alpha1,p2,p1);
    carea[c] += dca; darea[c] -= dca;
    if (j==ii) break;
    i=j; i1=j1;
    j1=perconcat(j1,vper[j][1]);
    j=vnbr[j][1];
    x2=x3; y2=y3;
  } while (TRUE);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     * C E L L P O P ( )    *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	cellpop(short i, k, REAL dptol)
 *
 *	Arguments:	i	= index of vertex adjacent to given cell
 *			k	= selects the k'th cell adjacent to `i'
 *			dptol	= patch the cell up sufficiently so that it
 *				is possible to adjust pressures by an amount
 *				up to `dptol' without breaking any arcs of the
 *				cell.
 *
 *	Return value:	none
 *
 *	Action:		This routine (along with its companion `bordpop()')
 *			is the fundamental solution to `arc breakage' (see
 *			the comments accompanying `barcangle()' for a brief
 *			discussion of this problem).
 *
 *			Whenever an arc breakage error condition is flagged,
 *			then this routine (or its companion `bordpop()') will
 *			ultimately be called to remedy the situation.
 *
 *			SPECIFICATION OF THE GIVEN CELL:
 *			The regions adjacent to vertex `i' are given in a
 *			standard order by the neighbour index `k = 0, 1 or 2'.  
 *			The value `k=0' always selects the adjacent Plateau
 *			border (an illegal choice for this routine) and `k=1'
 *			or `k=2' select the two adjacent cells in
 *			anti-clockwise order about vertex `i'.
 *
 *			ALGORITHM TO REMEDY BROKEN ARCS:
 *			The routine runs around all the arcs on the perimeter
 *			of the cell, `popping' any border arcs which are 
 *			broken.  But popping is only done in the direction that
 *			the arcs want to go.  For example if the cell is
 *			currently trying to increase its area (`darea[c]>0')
 *			then the arcs pop outwards, never inwards.  Once the
 *			arc popping has been done, the cell repairs the broken
 *			arcs by lowering its pressure.  The amount by which the
 *			cell pressure is lowered is chosen so as to allow room
 *			to manoeuvre for subsequent computation (e.g. the cell
 *			pressure is lowered by at least `dptol').
 *
 *****************************************************************************/
void cellpop(i, k, dptol)
short i, k;
REAL dptol;
{
  short i1, ii, ii1, j, j1, c, perconcat();
  REAL tx, ty, p1, p2, pb, pbmin, x1, y1, x2, y2, x3, y3, alpha, d, dsup, r,
       dpcrit;
  REAL barcangle(), carcangle(), barclen(), linlen(), fsign();
  boolean la, arcbrk, larc();
  void trans(), setlarc();
  if (k==2) { i=vnbr[i][0]; k=1; }
  c=cadj[i][k];
  ii=i; p1=cp[c];
  i1=0; j=vnbr[ii][2]; j1=vper[ii][2];
  tx=(REAL)PERX(j1); ty=(REAL)PERY(j1);
  x1=vx[ii]; y1=vy[ii]; x2=x1; y2=y1;
  dsup=0.0; arcbrk=FALSE; pbmin=0.0;
  do {
    x3=vx[j]+tx*boxwid; y3=vy[j]+ty*boxhgt;
    la=larc(i,2);
    /* Test for arcbreakage here */
    d=linlen(x2,y2,x3,y3); r=BRADIUS(p1,pb=bp[cadj[i][0]]);
    if (d > 2.0*r*cosminang) {
      /* Pop the arc! */
      if (darea[c]>0.0) {
        setlarc(&vper[i][2],la=TRUE); setlarc(&vper[j][1],la=TRUE);
      }
      else {
        setlarc(&vper[i][2],la=FALSE); setlarc(&vper[j][1],la=FALSE);
      }
      dsup= (dsup>d) ? dsup : d;
      pbmin= min(pbmin,pb);
      arcbrk=TRUE;
    }
    i=j; i1=j1;
    tx += PERX(vper[j][0]); ty += PERY(vper[j][0]);
    j=vnbr[j][0];
    x2=x3; y2=y3;
    x3=vx[j]+tx*boxwid; y3=vy[j]+ty*boxhgt;
    /* p2=cp[cadj[i][1]]; */
    if (j==ii) break;
    i=j; i1=j1;
    tx += PERX(vper[j][2]); ty += PERY(vper[j][2]);
    j=vnbr[j][2];
    x2=x3; y2=y3;
  } while (TRUE);
  /* Make sure that it is now safe to increment the cell pressure by */
  /* at least 'dptol', and safe to move all points on the boundary */
  /* at least the distance 'lengthdelta' */
  if (arcbrk) {
    dpcrit=1.0/BRADIUS(p1,pbmin)-2.0*cosminang/(dsup+lengthdelta);
    dpcrit= (dpcrit>dptol) ? dpcrit : dptol;
    cp[c] -= dpcrit;
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     * B O R D P O P ( )    *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	bordpop(short i, REAL dptol)
 *
 *	Arguments:	i	= index of vertex adjacent to Plateau border
 *			dptol	= patch the border up sufficiently so that it
 *				is possible to adjust pressures by an amount
 *				up to `dptol' without breaking any arcs of the
 *				Plateau border.
 *
 *	Return value:	none
 *
 *	Action:		This routine (along with its companion `cellpop()')
 *			is the fundamental solution to `arc breakage' (see
 *			the comments accompanying `barcangle()' for a brief
 *			discussion of this problem).
 *
 *			Whenever an arc breakage error condition is flagged,
 *			then this routine (or its companion `cellpop()') will
 *			ultimately be called to remedy the situation.
 *
 *			This routine runs around all the arcs on the perimeter
 *			of the Plateau border and, if it finds a broken arc,
 *			it `pops' the cell adjacent to this arc by calling
 *			`cellpop()' as a subroutine. (More details of what this
 *			entails can be got by reading the comments on the
 *			subroutine `cellpop()'.
 *
 *****************************************************************************/
void bordpop(i,dptol)
short i;
REAL dptol;
{
  short i1, ii, ii1, j, j1, b, perconcat();
  REAL r, d, p1, p2, x1, y1, x2, y2, x3, y3, fsign();
  void trans();
  b=cadj[i][0];
  ii=i; p1=bp[b];
  i1=0; j=vnbr[ii][1]; j1=vper[ii][1];
  x1=vx[ii]; y1=vy[ii]; x2=x1; y2=y1;
  do {
    trans(vx[j],vy[j],j1,&x3,&y3);
    p2=cp[cadj[i][2]];
    r=BRADIUS(p2,p1); d=linlen(x2,y2,x3,y3);
    if (d+lengthdelta > 2.0*r*cosminang) cellpop(i,2,dptol);
    if (j==ii) break;
    i=j; i1=j1;
    j1=perconcat(j1,vper[j][1]);
    j=vnbr[j][1];
    x2=x3; y2=y3;
  } while (TRUE);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  B S I D E S ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	short bsides(short b)
 *
 *	Arguments:	b	= index of Plateau border
 *
 *	Return value:	Number of sides of the Plateau border.
 *			(Note that this function is not particularly efficient
 *			since it has to search for a vertex which lies on
 *			this Plateau border.  It would probably make more
 *			sense to pass it a vertex index directly).
 *
 *****************************************************************************/
short bsides(b)
short b;
{
  short ii, i, nbs;
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    if (b==cadj[i][0]) break;
  }
  nbs=0; ii=i;
  do { i=vnbr[i][1]; nbs++;
  } while (i!=ii);
  return nbs;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  C S I D E S ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	short csides(short c)
 *
 *	Arguments:	c	= index of cell
 *
 *	Return value:	Number of sides of the cell.
 *			(Note that this function is not particularly efficient
 *			since it has to search for a vertex which lies on
 *			the cell.  It would probably make more sense to pass it
 *			a vertex index directly).
 *
 *****************************************************************************/
short csides(c)
short c;
{
  short i, ii, k, ncs;
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    if (c==cadj[i][k=1]) break;
    if (c==cadj[i][k=2]) break;
  }
  if (k==2) { i=vnbr[i][0]; k=1; }
  ncs=0; ii=i;
  do { i=vnbr[i][2]; i=vnbr[i][0]; ncs++;
  } while (i!=ii);
  return ncs;
}


#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/1000000000.0)

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *    R A N 3 ( )       *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL ran3(int idum)
 *
 *	Arguments:	idum	= dummy value when positive or, when negative,
 *				seeds the random number generator with |idum|.
 *
 *	Return value:	Random number.
 *
 *	Acknowledgements: Taken from the book `Numerical Recipes in C' by
 *			  Press, Flannery and Teukolsky, Cambridge.
 *
 *****************************************************************************/
REAL ran3(idum)
int idum;
{
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;

  if (idum < 0 || iff==0) {
    iff=1;
    mj=MSEED-( (idum<0) ? -idum : idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1; i<=54; i++) {
      ii=(21*i)%55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk<MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1; k<=4; k++)
      for (i=1; i<=55; i++) {
        ma[i] -= ma[1+(i+30)%55];
        if (ma[i]<MZ) ma[i] += MBIG;
      }
    inext=0;
    inextp=31;
    idum=1;
  }
  if (++inext==56) inext=1;
  if (++inextp==56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj<MZ) mj += MBIG;
  ma[inext]=mj;
  return (mj*FAC);
}
               
/**************************/
/* Error handling routine */
/**************************/

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     * P L E R R O R ( )    *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	plerror(char *s)
 *
 *	Arguments:	*s	= error message string
 *
 *	Return value:	none
 *
 *	Action:		Print an error message.
 *
 *****************************************************************************/
void plerror(s)
char *s;
{
  fprintf(stderr,"plateau: %s\n", s);
  fflush(stderr);
}
