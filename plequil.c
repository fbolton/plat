#include "include.h"
#ifdef IRIS
#include <device.h>
#endif

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     * P T R E L A X ( )    *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *		T O   B E   M O D I F I E D ! ! ! !
 *
 *	Subroutine:	ptrelax(short i, REAL *dx1, *dy1)
 *
 *	Arguments:	i	= vertex to be relaxed
 *			(*dx1, *dy1) = suggested displacement to be added to
 *				vertex `i' in order to move it towards an
 *				equilibrium position
 *
 *	Return value:	none
 *
 *	Action:		This subroutine carries out a fundamental step in
 *			relaxing a vertex so that it obeys the wet foam
 *			equilibrium conditions.  It `suggests' that the vertex
 *			`i' be displaced by the amount (*dx1, *dy1), in order
 *			come closer to an equilibrium configuration but it
 *			does not actually move the vertex (that is taken care
 *			of in the separate routine `vincrement()' called later
 *			on).
 *
 *			The magnitude of the suggested displacement is highly
 *			unreliable however it invariably tends to point in the
 *			right direction.  For this reason, the routine
 *			`vincrement()', which later actually moves the vertex,
 *			tends to limit the size of the suggested increment
 *			(*dx1, *dy1) while keeping the same direction.
 *
 *			ALGORITHM:
 *			This routine tries to estimate a displacement of the
 *			vertex `i' which would place it in a position where the
 *			two border arcs incident at `i' line up smoothly with
 *			the cell-cell arc incident at `i'.
 *
 *
 *****************************************************************************/
void ptrelax(i,dx1,dy1)
short i;
REAL *dx1, *dy1;
{
  short c1, c2, ii, j, k;
  REAL del, del2, x1, y1, x2, y2, x11, y11, x12, y12, p1, p2, pb, delp2,
       dplus, dminus, xc1, yc1, xc2, yc2, rx1, ry1, rx12, ry12, alpha1, alpha2,
       r12, rdot;
  REAL bcangle1(), bcangle2(), barcangle(), v[2], m[2][2], linlen();
  boolean la11, la12, larc(), arccentre(), arcbrk;
  void matinv2(), vnbrxy(), bordpop();

  del=lengthdelta; del2=del/2.0; delp2=pressuredelta/2.0;
  x1=vx[i]; y1=vy[i];
  vnbrxy(i,0,&x2,&y2); vnbrxy(i,1,&x11,&y11); vnbrxy(i,2,&x12,&y12);

  c1=cadj[i][2]; c2=cadj[i][1];
  p1=cp[c1]; p2=cp[c2]; pb=bp[cadj[i][0]];
  la11=larc(i,1); la12=larc(i,2);
  /* Now set up the vector angle discrepancies */
  v[0]= PI-bcangle1(x12,y12,x1,y1,x2,y2,p2,p1,pb,la12,&arcbrk);
  if (arcbrk) {
    bordpop(i,delp2); la11=larc(i,1); la12=larc(i,2);
    p1=cp[c1]; p2=cp[c2]; pb=bp[cadj[i][0]];
    v[0]= PI-bcangle1(x12,y12,x1,y1,x2,y2,p2,p1,pb,la12,&arcbrk);
  }
  v[1]= PI-bcangle2(x11,y11,x1,y1,x2,y2,p1,p2,pb,la11,&arcbrk);
  if (arcbrk) {
    bordpop(i,delp2); la11=larc(i,1); la12=larc(i,2);
    p1=cp[c1]; p2=cp[c2]; pb=bp[cadj[i][0]];
    v[1]= PI-bcangle2(x11,y11,x1,y1,x2,y2,p1,p2,pb,la11,&arcbrk);
  }
  if (v[1]>v[0]) {
    alpha1=barcangle(x1,y1,x11,y11,p1,pb,la11,&arcbrk);
    alpha2=barcangle(x1,y1,x12,y12,pb,p2,la12,&arcbrk);
    arccentre(x1,y1,x11,y11,alpha1,&xc1,&yc1);
    arccentre(x1,y1,x12,y12,alpha2,&xc2,&yc2);
    rx1=x1-xc1; ry1=y1-yc1;
    rx12=xc2-xc1; ry12=yc2-yc1; r12=sqrt(rx12*rx12+ry12*ry12);
    rx12 /= r12; ry12 /= r12;
    rdot=rx1*rx12+ry1*ry12;
    vx[i] += (rdot*rx12+xc1)-x1; vy[i] += (rdot*ry12+yc1)-y1;
    x1=vx[i]; y1=vy[i];
    v[0]= PI-bcangle1(x12,y12,x1,y1,x2,y2,p2,p1,pb,la12,&arcbrk);
    if (arcbrk) {
      bordpop(i,delp2); la11=larc(i,1); la12=larc(i,2);
      p1=cp[c1]; p2=cp[c2]; pb=bp[cadj[i][0]];
      v[0]= PI-bcangle1(x12,y12,x1,y1,x2,y2,p2,p1,pb,la12,&arcbrk);
    }
    v[1]= PI-bcangle2(x11,y11,x1,y1,x2,y2,p1,p2,pb,la11,&arcbrk);
    if (arcbrk) {
      bordpop(i,delp2); la11=larc(i,1); la12=larc(i,2);
      p1=cp[c1]; p2=cp[c2]; pb=bp[cadj[i][0]];
      v[1]= PI-bcangle2(x11,y11,x1,y1,x2,y2,p1,p2,pb,la11,&arcbrk);
    }
  }

  if (linlen(x1,y1,x2,y2)<minvvlen) {
    *dx1= *dy1=0.0; return;
  }

  dplus=bcangle1(x12,y12,x1+del2,y1,x2,y2,p2,p1,pb,la12,&arcbrk);
  if (arcbrk) {
    bordpop(i,delp2); la11=larc(i,1); la12=larc(i,2);
    p1=cp[c1]; p2=cp[c2]; pb=bp[cadj[i][0]];
    dplus=bcangle1(x12,y12,x1+del2,y1,x2,y2,p2,p1,pb,la12,&arcbrk);
  }
  dminus=bcangle1(x12,y12,x1-del2,y1,x2,y2,p2,p1,pb,la12,&arcbrk);
  if (arcbrk) {
    bordpop(i,delp2); la11=larc(i,1); la12=larc(i,2);
    p1=cp[c1]; p2=cp[c2]; pb=bp[cadj[i][0]];
    dminus=bcangle1(x12,y12,x1-del2,y1,x2,y2,p2,p1,pb,la12,&arcbrk);
  }
  m[0][0]=(dplus-dminus)/del;
  dplus=bcangle2(x11,y11,x1+del2,y1,x2,y2,p1,p2,pb,la11,&arcbrk);
  if (arcbrk) {
    bordpop(i,delp2); la11=larc(i,1); la12=larc(i,2);
    p1=cp[c1]; p2=cp[c2]; pb=bp[cadj[i][0]];
    dplus=bcangle2(x11,y11,x1+del2,y1,x2,y2,p1,p2,pb,la11,&arcbrk);
  }
  dminus=bcangle2(x11,y11,x1-del2,y1,x2,y2,p1,p2,pb,la11,&arcbrk);
  if (arcbrk) {
    bordpop(i,delp2); la11=larc(i,1); la12=larc(i,2);
    p1=cp[c1]; p2=cp[c2]; pb=bp[cadj[i][0]];
    dminus=bcangle2(x11,y11,x1-del2,y1,x2,y2,p1,p2,pb,la11,&arcbrk);
  }
  m[1][0]=(dplus-dminus)/del;
  dplus=bcangle1(x12,y12,x1,y1+del2,x2,y2,p2,p1,pb,la12,&arcbrk);
  if (arcbrk) {
    bordpop(i,delp2); la11=larc(i,1); la12=larc(i,2);
    p1=cp[c1]; p2=cp[c2]; pb=bp[cadj[i][0]];
    dplus=bcangle1(x12,y12,x1,y1+del2,x2,y2,p2,p1,pb,la12,&arcbrk);
  }
  dminus=bcangle1(x12,y12,x1,y1-del2,x2,y2,p2,p1,pb,la12,&arcbrk);
  if (arcbrk) {
    bordpop(i,delp2); la11=larc(i,1); la12=larc(i,2);
    p1=cp[c1]; p2=cp[c2]; pb=bp[cadj[i][0]];
    dminus=bcangle1(x12,y12,x1,y1-del2,x2,y2,p2,p1,pb,la12,&arcbrk);
  }
  m[0][1]=(dplus-dminus)/del;
  dplus=bcangle2(x11,y11,x1,y1+del2,x2,y2,p1,p2,pb,la11,&arcbrk);
  if (arcbrk) {
    bordpop(i,delp2); la11=larc(i,1); la12=larc(i,2);
    p1=cp[c1]; p2=cp[c2]; pb=bp[cadj[i][0]];
    dplus=bcangle2(x11,y11,x1,y1+del2,x2,y2,p1,p2,pb,la11,&arcbrk);
  }
  dminus=bcangle2(x11,y11,x1,y1-del2,x2,y2,p1,p2,pb,la11,&arcbrk);
  if (arcbrk) {
    bordpop(i,delp2); la11=larc(i,1); la12=larc(i,2);
    p1=cp[c1]; p2=cp[c2]; pb=bp[cadj[i][0]];
    dminus=bcangle2(x11,y11,x1,y1-del2,x2,y2,p1,p2,pb,la11,&arcbrk);
  }
  m[1][1]=(dplus-dminus)/del;
  matinv2(m);
  *dx1=m[0][0]*v[0]+m[0][1]*v[1];
  *dy1=m[1][0]*v[0]+m[1][1]*v[1];
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     * M A T I N V ( )      *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	matinv2(REAL m[2][2])
 *
 *	Arguments:	m[2][2]	= a 2x2 real matrix
 *
 *	Return value:	none
 *
 *	Action:		Inverts the given 2x2 real matrix (or returns an error
 *			if the matrix is uninvertible).
 *			Used by the subroutine `ptrelax()'.
 *
 *****************************************************************************/
void matinv2(m)
REAL m[2][2];
{
  REAL x;
  void plerror();
  x=m[0][0]; m[0][0]=m[1][1]; m[1][1]=x;
  m[0][1]= -m[0][1]; m[1][0]= -m[1][0];
  x=m[0][0]*m[1][1]-m[0][1]*m[1][0];
  if (x==0.0) {
    plerror("uninvertible 2x2 matrix");
    m[0][0]=0.0; m[0][1]=0.0; m[1][0]=0.0; m[1][1]=0.0;
  }
  else {
    m[0][0] /= x; m[0][1] /= x; m[1][0] /= x; m[1][1] /= x;
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *   E Q U I L ( )      *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL equil(boolean adjustflag)
 *
 *	Arguments:	adjustflag	= when set to TRUE, an `adjustment' is
 *				made to vertices instead of a full
 *				equilibration (usually passed as FALSE).
 *
 *	Return value:	Returns `sup' which is the largest increment given to
 *			any vertex during this equilibration step (thus it can
 *			be used to estimate how close the network has got to
 *			equilibrium).
 *
 *	Action:		This subroutine does one iteration of a froth
 *			equilibration.  Here equilibrium is defined by two
 *			conditions:
 *			(i) The area of each cell `k' must match up with the
 *			target area `carea[k]+darea[k]'
 *			(ii) The slopes of all the arcs meeting at a point
 *			must be equal.
 *
 *			The program also has an optional effect, selected by
 *			`adjustflag = TRUE', whereby the vertices will all be
 *			shifted by the amounts given in (dvx[], dvy[]) and
 *			some housekeeping tasks will be done.  No equilibration
 *			is done in this mode.
 *
 *			Usually the subroutine is called as `equil(0)' and the
 *			equilibration algorithms are also carried out.  The
 *			algorithm carries out the following steps:
 *			(i) Area equilibration:  Firstly the areas of the
 *			cells are equilibrated so that they equal
 *			`carea[i]+darea[i]' to a good degree of accuracy.
 *			This is done primarily by adjusting the pressures
 *			`cp[i]'.  However, when a froth is very dry it is also
 *			necessary to adjust the pressures of the Plateau
 *			borders `bp[j]' -- the case of a dry froth is
 *			signalled by the flag `foamlike'.
 *			(ii) Fix the new Plateau border areas to be equal to
 *			their actual values
 *			(iii) Vertex equilibration: A *single* increment is
 *			worked out as a suggested change in position of the
 *			vertices (note that unlike `area equilibration' this is
 *			not iterated within this subroutine).
 *			(iv) Edge loss detection:  It is worked out if the
 *			suggested movements of the vertices will give rise
 *			to any edges being lost -- these edges are recorded in
 *			`iel[]'.
 *			(v) Vertex incrementation: Then increments are given
 *			to the vertices by calling the subroutine
 *			`vincrement()'.  The reason for being circumspect
 *			about this is that you can break the network by moving
 *			vertices.
 *			(vi) Border pinching:  Now we deal with the first kind
 *			of topological change.  If a Plateau border has more
 *			than three sides then two of its arcs can
 *			`pinch' together, thereby splitting off a new Plateau
 *			border.  The subroutine `bpinch()' is called here both
 *			to test for and implement this change.
 *			(vii) Edge loss:  The list of edges which were
 *			recorded in `iel[]' are now deleted.  The nett effect
 *			of this is that Plateau borders will coalesce.
 *			(viii) Cell area update:  The new geometrical areas
 *			of all the cells are recalculated and stored in
 *			`carea[]' (note that the target cell areas = 
 *			`carea[i]+darea[i]' remain the same).
 *
 *
 *****************************************************************************/
REAL equil(adjustflag)
boolean adjustflag;
{
/* This drives the equilibration and topological routines, doing
a single pass of the whole network.  The value returned is a figure to
indicate how converged the network is.
   The convergence error used is the supremum of the increments to the
vertices. */
  short i, i1, i2, k, j, j1, k1, c1, c2, b, ii, icount;
  REAL x1,y1,x11,y11,x12,y12,x2,y2,x21,y21,x22,y22,p1,p2,da1,da2,dp1,dp2,sup;
  REAL prate, vrate, vvrate, maxdcp, maxcp, dcp1, delp, delp2, dadp, maxdafrac,
       bmaxdafrac, daconverge, omaxdafrac=0.0, f, ca, cl;
  REAL a1, a2, a3, dx, dy, lim, d, d1;
  boolean la11, la12, la21, la22, arcbrk, stuck;
  boolean found[MBORD];
  void vnbrxy(), trans(), erelax(), elose(), bpinch(), normalizep(),
       normalizeba();
  void plerror(), ptrelax(), careaperim(), dvlimit(), cellpop(), bordpop();
  short perconcat();
  REAL cellarea(), bordarea(), fsign(), linlen();
  boolean larc(), elost(), vincrement();
  void foamplot();
  /*
   * Initialise:
   *
   * nel = number of edges to lost due to topological changes
   * sup = size of largest suggested increment to a vertex
   */
  nel=0; sup=0.0;
  /*
   * When the flag `adjustflag == TRUE' then the routine does *not* do
   * a full equilibration.  Instead it just increments the vertices by
   * the amounts given in (dvx[], dvy[]) and does some housekeeping to
   * make sure that the topology/geometry is not broken by these changes.
   *
   * This option is useful, for example, when you want to subject the whole
   * network to Hencky strain -- the vertices of the network can be subjected
   * to a linear transformation before the equilibration is begun.
   */
  if (adjustflag) goto ADJUST;
  /*
   * Shorthand for a `small' amount of pressure delta p...
   * used for calculating numerical derivatives w.r.t pressure etc.
   */
  delp=pressuredelta; delp2=delp/2.0;
  icount=0; maxdafrac=0.0; daconverge=0.0; prate=spdamp;
  /*
   * Adjust the total area of the network so that it fits within the given
   * periodic box.  This is just a housekeeping task in case this total
   * area has drifted due to cumulative numerical error.
   */
  normalizeba();
  /*
   * `do {} while' Loop for AREA EQUILIBRATION...
   *
   * The first step of equilibration is to adjust the pressures of all
   * the cells and Plateau borders until all of the cell areas are
   * equal to their target areas.  This condition is imposed to a fairly high
   * degree of accuracy, as experience shows it needs to be.
   */
  do {
    /*
     * This loop can get stuck if the suggested pressure increments
     * get too small.  Use this flag to avoid getting trapped within this
     * loop forever ! ...
     */
    stuck=TRUE;
    bmaxdafrac=0.0;
    /*
     * A froth is defined to be `foamlike' (actually, what is meant is `like
     * a dry foam') if the Plateau borders are very small.  The exact
     * definition of `foamlike' is given in the routine `setconstants()'.
     *
     * When the Plateau borders are relatively small, the simulation behaves
     * somewhat differently -- in particular here it is necessary to
     * individually equilibrate Plateau borders (which is something we don't
     * bother with when the borders are large).
     */
    if (foamlike) {
      for (i=0; i<onb; i++) found[i]=FALSE;
      /*
       * Loop around all Plateau borders and do:
       * Numerically calculate the derivative
       * dadp = d{area of Plateau border} / d{pressure of Plateau border}
       * ...and use this to estimate an increment to border pressure via
       * suggested dp = ({ideal border area} - {actual border area}) / dadp
       *
       * Note however, that {ideal border area} is not absolutely fixed (that
       * would be physically unrealistic) but is reset to its actual value
       * later on in this subroutine.

       * Also calculate `bmaxdafrac', the maximum fractional change in area of
       * a Plateau border, as a convergence criterion for area equilibration.
       */
      for (ii=0; ii<nv; ii++) {
        i=vlist[ii];
        if (!found[b=cadj[i][0]]) {
          found[b]=TRUE;
          bp[b] += delp2;
          a3=bordarea(i,&arcbrk);
          if (arcbrk) {
            bordpop(i,delp2); a3=bordarea(i,&arcbrk);
          }
          bp[b] -= delp;
          a1=bordarea(i,&arcbrk);
          if (arcbrk) {
            bordpop(i,delp2); a1=bordarea(i,&arcbrk);
          }
          bp[b] += delp2;
          a2=bordarea(i,&arcbrk);
          if (arcbrk) {
            bordpop(i,delp2); a2=bordarea(i,&arcbrk);
          }
          dadp=(a3-a1)/delp;
          dp1=0.5*sbpdamp*(barea[b]-a2)/dadp;
          bmaxdafrac=
                 (bmaxdafrac>(f=fabs((a2-barea[b])/barea[b])))?bmaxdafrac : f;
          if (fabs(dp1)>delp) {
            stuck=FALSE;
            /* Restrict increments to 0.1 of 'bpav' */
            if (fabs(dp1/bpav)>0.1) dp1= fsign(dp1)*0.1*fabs(bpav);
            bp[b] += dp1;
          }
        }
      }
    }
    /*
     * Clear array `found[]' *and* array `dcp[]' which will store
     * suggested changes in cell pressure (note that the suggested
     * changes to a cell's pressure will come not just from the cell itself
     * but also from its immediate neighbours).
     */
    for (i=0; i<onc; i++) { found[i]=FALSE; dcp[i]=0.0; }
    /*
     * Loop around all the cells:
     *
     * Calculate the numerical derivative of each cell
     * dadp = d{area of cell} / d{pressure of cell}
     * ...and then estimate the amount by which to change pressure
     * suggested dcp = 0.5*({target area} - {actual area})/dadp
     * The factor of `0.5' comes in because the cell itself has its
     * pressure *increased* by half the amount, while its neighbours
     * have their pressure *decreased* so as to account for the rest
     * (roughly).
     */
    for (ii=0; ii<nv; ii++) {
      i=vlist[ii];
      for (k=1; k<3; k++) {
        if (!found[c1=cadj[i][k]]) {
          found[c1]=TRUE;
          a2=cellarea(i,k,&arcbrk);
          if (arcbrk) {
            cellpop(i,k,delp2); a2=cellarea(i,k,&arcbrk);
          }
          darea[c1]=carea[c1]+darea[c1]-a2;
          carea[c1]=a2;
          /* a2=carea[c1]; */
          cp[c1] += delp2;
          a3=cellarea(i,k,&arcbrk);
          if (arcbrk) {
            cellpop(i,k,delp2); a3=cellarea(i,k,&arcbrk);
          }
          cp[c1] -= delp;
          a1=cellarea(i,k,&arcbrk);
          if (arcbrk) {
            cellpop(i,k,delp2); a1=cellarea(i,k,&arcbrk);
          }
          cp[c1] += delp2;
          dadp=(a3-a1)/delp;
          if (dadp<0.0) {
            /* plerror("negative area derivative in routine equil"); */
            /* this is actually correct in some cases where a pb arc */
            /* has popped. */
          }
          dp1=0.5*darea[c1]/dadp;
          /* the factor of 0.5 comes in here because it is assumed that, */
          /* on average, its nbrs provide the other half of the required */
          /* change in pressure difference. */
          if (/* (-cp[c1]/bpav)<0.1 && */ fabs(dp1)>delp) {
            stuck=FALSE;
            dcp[c1] += dp1;
            /* Now change the pressure of its nbrs by a complementary amount */
            dp1= -dp1/((REAL) ncsides[c1]);
            if (k==2) { i=vnbr[i][0]; k=1; }
            j=i;
          }	
        }
      }
    }
    /* Now to determine what fraction of the calculated increments ought */
    /* to be applied */
    /*
     * First get the largest cell pressure increment...
     */
    maxdcp=0.0;
    for (ii=0; ii<nc; ii++) {
      i=clist[ii];
      maxdcp= (maxdcp > (dcp1=fabs(dcp[i]))) ? maxdcp : dcp1;
    }
    /*
     * ... then restrict the largest pressure change to 0.01 of 'bpav'
     */
    if (fabs(maxdcp/bpav)<0.01) prate=spdamp;
    else {
      prate=spdamp*fabs(0.01*bpav/maxdcp);
    }
    /*
     * Now decrease the pressure increments by 'prate' and add them
     * to the cell pressures `cp[]'
     */
    maxdcp *= prate;
    for (ii=0; ii<nc; ii++) {
      i=clist[ii];
      dcp[i] *= prate;
      cp[i] += dcp[i];
    }
    /*
     * Now correct the areas of the cells...
     *
     * ...and calculate `maxdafrac', the maximum fractional change in area,
     * as a convergence criterion for the area equilibration.
     */
    omaxdafrac=maxdafrac; maxdafrac=0.0; maxcp=0.0;
    for (i=0; i<onc; i++) found[i]=FALSE;
    for (ii=0; ii<nv; ii++) {
      i=vlist[ii];
      for (k=1; k<3; k++) {
        if (!found[c1=cadj[i][k]]) {
          found[c1]=TRUE;
          da1=cellarea(i,k,&arcbrk);
          if (arcbrk) { cellpop(i,k,delp2); da1=cellarea(i,k,&arcbrk); }
          da1 -= carea[c1];
          carea[c1] += da1; darea[c1] -= da1;
          maxdafrac= (maxdafrac>(f=fabs(darea[c1]/carea[c1]))) ? maxdafrac : f;
          maxcp= (maxcp>(p1=cp[c1])) ? maxcp : p1;
        }
      }
    }
    daconverge= (omaxdafrac==0.0) ? 0.0 : maxdafrac/omaxdafrac;
    icount++;
    // foamplot(0);
  } while (((maxdafrac>areasup) || (bmaxdafrac>bareasup)) && !stuck);
  /*
   * ...finished AREA EQUILIBRATION.
   *
   * Now shift all pressures in the system so that the average cell area
   * works out to be zero (just a handy convention).
   */
  normalizep();
  /*
   * Fix the new border areas at the present values
   * (as we promised we would earlier on)
   */
  for (i=0; i<onb; i++) found[i]=FALSE;
  for (ii=0; ii<nb; ii++) {
    i=vlist[ii];
    if (!found[b=cadj[i][0]]) {
      found[b]=TRUE;
      da1=bordarea(i,&arcbrk);
      if (arcbrk) { bordpop(i,delp2); da1=bordarea(i,&arcbrk); }
      barea[b]=da1;
    }
  }
  /*
   * Now do VERTEX EQUILIBRATION...
   *
   * The key task is performed by `ptrelax()' which suggests
   * the increments (dvx[], dvy[]) which should be made to
   * relax a vertex.
   */
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    ptrelax(i,&dvx[i],&dvy[i]);
  }
ADJUST : ;
  if (adjustflag) { vrate=1.0; vvrate=1.0; }
  else if (foamlike) { vrate=0.5*svdamp; vvrate=0.5*svvdamp; }
  else { vrate = vvrate = 0.5*svvdamp; }
  /* the factor of 0.5 anticipates the fact that the increments */
  /* tend to overshoot */
  /*
   * Calculate the convergence criterion `sup'...
   * (which is the largest suggested increment to a vertex)
   */
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    sup= (sup>fabs(dvx[i])) ? sup : fabs(dvx[i]);
    sup= (sup>fabs(dvy[i])) ? sup : fabs(dvy[i]);
  }
  if (adjustflag) goto ADJUST1;
  for (i=0; i<onb; i++) found[i]=FALSE;
  /*
   * Turn the suggested increments (dvx[], dvy[]) into reasonable
   * sized increments.  As was noted in the routine `ptrelax()', the
   * increments which are returned are often unreasonably large
   * (e.g. when the vertex is far from equilibrium) although they tend
   * to point in the right direction.  Therefore, we limit the size of
   * the increments here...
   *
   * In particular, `vvrate' is the usual rate of damping of a vertex
   * increment.  Additionally, an absolute upper limit to the increment
   * size is imposed by the subroutine `dvlimit()'.
   * 
   * Note that three-sided Plateau borders get special treatment.  There
   * are two rates of damping used -- `vrate', gentle damping, for the
   * motion of the centre of mass of the 3-border; and `vvrate', strong
   * damping, used for the relative motion of the vertices within the
   * 3-border.  This allows networks with small 3-sided Plateau borders
   * to equilibrate more rapidly.
   */
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    if (!found[b=cadj[i][0]]) {
      found[b]=TRUE;
      if (nbsides[b]>3) {
        dvx[i] *= vvrate; dvy[i] *= vvrate;
        dvlimit(maxdvv,i,dvx,dvy);
      }
      else {  /* 3-sided borders get special treatment! */
        /* Essentially, the average motion of the 3-border is damped */
        /* by 'vrate',  while the distortion of the 3-border is */
        /* damped by 'vvrate'. */
        i1=vnbr[i][1]; i2=vnbr[i][2];
        dx=(dvx[i]+dvx[i1]+dvx[i2])/3.0;
        dy=(dvy[i]+dvy[i1]+dvy[i2])/3.0;
        dvx[i] -= dx; dvx[i1] -= dx; dvx[i2] -= dx;
        dvy[i] -= dy; dvy[i1] -= dy; dvy[i2] -= dy;
        dx *= vrate; dy *= vrate;
        lim=maxdv;
        if (fabs(dx)>fabs(dy)) {
          if (fabs(dx)>lim) {
            dx = fsign(dx)*lim;
            dy *= lim/fabs(dx);
          }
        }
        else {
          if (fabs(dy)>lim) {
            dx *= lim/fabs(dy);
            dy = fsign(dy)*lim;
          }
        }
        dvx[i] *= vrate; dvx[i1] *= vrate; dvx[i2] *= vrate;
        dvy[i] *= vrate; dvy[i1] *= vrate; dvy[i2] *= vrate;
        dvlimit(maxdvv,i,dvx,dvy);
        dvx[i] += dx; dvx[i1] += dx; dvx[i2] += dx;
        dvy[i] += dy; dvy[i1] += dy; dvy[i2] += dy;
      }
    }
  }
#ifndef BACKGROUND
  if (sup>0.2) {
    /* printf("warning: very large increment occurred\n"); */
  }
#endif
ADJUST1 : ;
  /*
   * Now calculate a suitable increment for the vertices of a one-sided cell.
   * They are moved more or less in parallel to the their neighbouring
   * vertices (ensuring that the one-sided cell remains attatched to its
   * neighbour).
   */
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    if (ncsides[cadj[i][1]]==1) {
      j=vnbr[i][0]; i1=vnbr[i][1]; j1=vnbr[j][2];
      x1=vx[i]; y1=vy[i]; vnbrxy(i,1,&x11,&y11);
      vnbrxy(i,0,&x2,&y2);
      k=perconcat(vper[i][0],vper[j][2]);
      trans(vx[j1],vy[j1],k,&x22,&y22);
      dx=dvx[j1]-dvx[i1]; dy=dvy[j1]-dvy[i1];
      dvx[i]=dvx[j]=dvx[i1]; dvy[i]=dvy[j]=dvy[i1];
      d=linlen(x11,y11,x22,y22); d1=linlen(x11,y11,x1,y1);
      f=d1/d;
      dvx[i] += f*dx; dvy[i] += f*dy;
      d1=linlen(x11,y11,x2,y2);
      f=d1/d;
      dvx[j] += f*dx; dvy[j] += f*dy;
    }
  }
  /*
   * Now test for lost edges
   */
  nel=0;
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    j=vnbr[i][0];
    if (j>i) {
      /* avoids double counting */
      x1=vx[i]; y1=vy[i]; vnbrxy(i,0,&x2,&y2);
      if (elost(x2-x1,y2-y1,dvx[j]-dvx[i],dvy[j]-dvy[i])) {
        /*
	 * Add it to the list of cell edges to be deleted...
	 */
        iel[nel]=i; nel++;
      }
    }
  }
  /*
   * At last!  Add the increments to the vertices...
   */
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    /* takes care of some knotty eventualities, even calling itself. */
    vincrement(i,dvx[i],dvy[i]);
  }
  /*
   * Loop over all Plateau borders, calling routine `bpinch()' for each
   * one.  Routine `bpinch()' tests for the second kind of topological
   * change, where two arcs of a Plateau border `pinch' together, thus
   * splitting off a piece of the Plateau border.
   */
  if (!notopol) {
    for (b=0; b<onb; b++) found[b]=FALSE;
    for (ii=0; ii<nv; ii++) {
      i=vlist[ii];
      b=cadj[i][0];
      if (!found[b]) {
        found[b]=TRUE;
        if (nbsides[b]>3) bpinch(i);
      }
    }
  }
  /*
   * Now implement the lost edges.
   * NB: a maximum of one edge adjacent to a given cell is lost,
   * since this aids stability
   */
  if (!notopol) {
    for (i=0; i<onc; i++) found[i]=FALSE;
    for (ii=0; ii<nel; ii++) {
      i=iel[ii];
      if (!found[cadj[i][1]] && !found[cadj[i][2]]) {
        found[cadj[i][1]]=TRUE; found[cadj[i][2]]=TRUE;
        elose(i);
      }
    }
  }
  /*
   * Now recalculate the cell areas and calculate the network energy
   * by adding up the total perimeter length of the froth...
   */
  maxdafrac=0.0; netenergy=0.0;
  for (i=0; i<onc; i++) found[i]=FALSE;
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    for (k=1; k<3; k++) {
      if (!found[c1=cadj[i][k]]) {
        found[c1]=TRUE;
        careaperim(i,k,TRUE,TRUE,&ca,&cl,&arcbrk);
        if (arcbrk) {
          cellpop(i,k,delp2); careaperim(i,k,TRUE,TRUE,&ca,&cl,&arcbrk);
        }
        netenergy += cl;
        da1=ca-carea[c1];
        carea[c1] += da1; darea[c1] -= da1;
        maxdafrac= (maxdafrac>(f=fabs(darea[c1]/carea[c1]))) ? maxdafrac : f;
      }
    }
  }
  return sup;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     * D V L I M I T ( )    *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	dvlimit(REAL lim, short i, REAL dvx[], dvy[])
 *
 *	Arguments:	lim	= size of maximum allowed vertex increment
 *			i	= index of vertex lying on the particular
 *				Plateau border in question
 *			(dvx[], dvy[]) = array of vertex increments
 *
 *	Return value:	none
 *
 *	Action:		Limits the maximum size of a vertex increment to `lim'.
 *			In particular, it takes the Plateau border selected
 *			via vertex index `i'; runs around the vertices
 *			lying on the Plateau border and finds the size of the
 *			largest of all the increments; it then proportionally
 *			reduces the increments to all of the vertices so
 *			that the largest becomes equal to `lim'.
 *
 *			The reason for treating a whole Plateau border all in
 *			one go, instead of individually limiting increments
 *			to vertices, is that we want to avoid dramatically
 *			distorting a Plateau border.
 *
 *****************************************************************************/
void dvlimit(lim,i,dvx,dvy)
short i;
REAL lim, dvx[], dvy[];
{
  short j, je;
  REAL sup, m, f, fsign();
  /* First of all find the largest increment: */
  sup=0.0;
  je=j=i;
  do {
    m=max(fabs(dvx[j]),fabs(dvy[j]));
    sup=max(m,sup);
    j=vnbr[j][1];
  } while (j != je);
  /* Then impose 'lim' as the largest increment, shrinking all increments */
  /* in proportion */
  if (sup>lim) {
    f=lim/sup;
    je=j=i;
    do {
      dvx[j] *= f; dvy[j] *= f;
      j=vnbr[j][1];
    } while (j != je);
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *   E L O S T ( )      *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	elost(REAL rx, ry, dx, dy)
 *
 *	Arguments:	(rx, ry) = the vector joining the endpoints of edge
 *			(dx, dy) = proposed total increment to edge vector
 *
 *	Return value:	TRUE if the edge is to be lost, else FALSE
 *
 *	Action:		Test whether the given edge disappears once the
 *			proposed vertex increments are made.  If the component
 *			of the increment vector lying along the edge vector is
 *			both negative and bigger than the edge vector, then
 *			the edge is lost.  Alternatively if the edge has
 *			simply got very small (less than `minvvlen') then it
 *			is also lost.
 *
 *****************************************************************************/
boolean elost(rx,ry,dx,dy)
REAL rx,ry,dx,dy;
{
  REAL rr;
  REAL minvv2=0.0;
  minvv2=minvvlen*minvvlen;
  rr=rx*rx+ry*ry;
  return (rr< -(rx*dx+ry*dy) || rr<minvv2);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *   E L O S E ( )      *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	elose(short k)
 *
 *	Arguments:	k	= index of one of the vertices on the edge
 *
 *	Return value:	none
 *
 *	Action:		First of all the two merging borders are moved
 *			together so as to close the small gap along the
 *			disappearing edge.  The algorithm then carries out the
 *			necessary bookkeeping to make the edge disappear.
 *
 *			There are two special cases which arise.  At the outset
 *			the routine checks if, by removing this edge, a
 *			percolating Plateau border would be created.  The
 *			program stops if this is the case.  The other special
 *			case is when either of the cells adjacent to the edge
 *			are single-sided.  In that case the edge is left alone
 *			since we do not want to cut off the cell completely
 *			from the network (which would oblige us to deal with
 *			a new object, the isolated bubble).
 *
 *****************************************************************************/
void elose(k)
short k;
{
  short i, ii, j, jj, kk, i1, j1, b, b1, c1, c2, ns, iseed, iseed1, iseed2;
  REAL dx, dy, x1, y1, x2, y2, p, pb, alpha, f, d, r;
  boolean la, la1, arcbrk;
  char s[80];
  void cforgetindex(), trans(), setlarc(), vnbrxy(), plerror(), bforgetindex(),
       vforgetindex(), setconstants(), bublistindex(), cellpop(), bordpop(),
       clusterplot(), foamout(), infoclose();
  short perconcat(), bgetindex();
  REAL barcangle(), linlen(), bordarea();
  boolean vincrement(), larc(), perculates();
  /*
   * The first special case:
   *
   * One-sided cells do not have their last edge removed,
   * and this is because we do not want to deal with isolated
   * bubbles.  So just return...
   */
  if (ncsides[cadj[k][1]]==1 || ncsides[cadj[k][2]]==1) return;
  /*
   * Set up the basic objects we are dealing with, for example
   * (x1, y1) and (x2, y2) are the endpoints of the edge to be removed
   */
  iseed1= -1; iseed2= -1;
  kk=vnbr[k][0];
  x1=vx[k]; y1=vy[k]; vnbrxy(k,0,&x2,&y2);
  b=cadj[k][0]; b1=cadj[kk][0];
  c1=cadj[k][1]; c2=cadj[k][2];
  pb=0.5*(bp[b]+bp[b1]);
  /*
   * This is the second special case:
   *
   * If the two Plateau borders which are about to merge are, 
   * in fact, one and the same Plateau border it indicates that
   * the froth is about to be broken apart by this percolating
   * Plateau border.
   */
  if (b==b1) { /* implies that the foam is breaking off a piece */
    if (perculates(k)) {
      sprintf(s,"perculation occured next to vertex %d",k); plerror(s);
      fflush(stdout);
      clusterplot(k);
      foamout("fmperc.pb");
#ifdef IRIS
      sleep(30);
#endif
      infoclose();
      exit(0);
    }
    else return;
  }
  /*
   * The usual case:
   *
   * Start by moving the two Plateau borders up alongside
   * each other so that we lose the edge gently.  Note that
   * the *entire* Plateau borders are moved so that the endpoints
   * lie on top of each other.
   */
  i=k; j=kk;
  dx=0.5*(x2-x1)-dvx[i]; dy=0.5*(y2-y1)-dvy[i];
  i1=i;
  do {
    dvx[i] += dx; dvy[i] += dy;
    i=vnbr[i][1];
  } while (i != i1);
  dx=0.5*(x1-x2)-dvx[j]; dy=0.5*(y1-y2)-dvy[j];
  j1=j;
  do {
    dvx[j] += dx; dvy[j] += dy;
    j=vnbr[j][1];
  } while (j != j1);
  /*
   * Deal with topological changes on the c1 side.
   *
   * This involves some straightforward bookkeeping to
   * join two Plateau border arcs together.
   */
  if (vnbr[k][2]==kk) {
    /* Do nothing! */ ;
    /*
     * This is where an isolated bubble would be created
     * if we were interested in making them!
     */
  }
  else {
    iseed=iseed1=i=vnbr[k][2]; ii=vnbr[kk][1];
    la=larc(i,1); la1=larc(kk,1);
    p=cp[c1];
    x1=vx[i]; y1=vy[i];
    i1=vper[i][1]; trans(vx[k],vy[k],i1,&x2,&y2);
    alpha=barcangle(x1,y1,x2,y2,p,pb,la,&arcbrk);
    if (arcbrk) {
      cellpop(k,1,pressuredelta/2.0);
      alpha=barcangle(x1,y1,x2,y2,p,pb,la,&arcbrk);
    }
    x1=x2; y1=y2;
    i1=perconcat(i1,vper[k][0]); i1=perconcat(i1,vper[kk][1]);
    trans(vx[ii],vy[ii],i1,&x2,&y2);
    alpha += barcangle(x1,y1,x2,y2,p,pb,la1,&arcbrk);
    if (arcbrk) {
      cellpop(k,1,pressuredelta/2.0);
      alpha += barcangle(x1,y1,x2,y2,p,pb,la1,&arcbrk);
    }
    /* Now do updates */
    vnbr[i][1]=ii; vnbr[ii][2]=i;
    vper[i][1]=i1;
    vper[ii][2]=PERFN(-PERX(i1),-PERY(i1));
    setlarc(&vper[i][1],(alpha>PI));
    setlarc(&vper[ii][2],(alpha>PI));
    x1=vx[i]; y1=vy[i];
    vnbrxy(i,1,&x2,&y2);
    d=linlen(x1,y1,x2,y2); r=BRADIUS(p,pb);
    if (2.0*r*cosminang<d) cellpop(i,2,pressuredelta/2.0);
  }
  /*
   * Deal with topological changes on the c2 side.
   *
   * This involves some straightforward bookkeeping to
   * join two Plateau border arcs together.
   */
  if (vnbr[k][1]==kk) {
    /* Do nothing! */ ;
    /*
     * This is where an isolated bubble would be created
     * if we were interested in making them!
     */
  }
  else {
    iseed=iseed2=j=vnbr[kk][2]; jj=vnbr[k][1];
    la=larc(j,1); la1=larc(k,1);
    p=cp[c2];
    x1=vx[j]; y1=vy[j];
    j1=vper[j][1]; trans(vx[kk],vy[kk], j1,&x2,&y2);
    alpha=barcangle(x1,y1,x2,y2,p,pb,la,&arcbrk);
    if (arcbrk) {
      cellpop(k,2,pressuredelta/2.0);
      alpha=barcangle(x1,y1,x2,y2,p,pb,la,&arcbrk);
    }
    x1=x2; y1=y2;
    j1=perconcat(j1,vper[kk][0]); j1=perconcat(j1,vper[k][1]);
    trans(vx[jj],vy[jj],j1,&x2,&y2);
    alpha += barcangle(x1,y1,x2,y2,p,pb,la1,&arcbrk);
    if (arcbrk) {
      cellpop(k,2,pressuredelta/2.0);
      alpha += barcangle(x1,y1,x2,y2,p,pb,la1,&arcbrk);
    }
    /* Now do updates */
    vnbr[j][1]=jj; vnbr[jj][2]=j;
    vper[j][1]=j1;
    vper[jj][2]=PERFN(-PERX(j1),-PERY(j1));
    setlarc(&vper[j][1],(alpha>PI));
    setlarc(&vper[jj][2],(alpha>PI));
    x1=vx[j]; y1=vy[j];
    vnbrxy(j,1,&x2,&y2);
    d=linlen(x1,y1,x2,y2); r=BRADIUS(p,pb);
    if (d>2.0*r*cosminang) cellpop(j,2,pressuredelta/2.0);
  }
  /*
   * Now do common tasks - ie eliminate vertices 'k', 'kk' etc.
   */
  vforgetindex(k); vforgetindex(kk);
  (ncsides[c1])--; (ncsides[c2])--;
  if (b!=b1) {
    /* Choose the smaller border index as the index of the new
    co-alesced border */
    if (b>b1) { k=b; b=b1; b1=k; }
    nbsides[b] += nbsides[b1]-2;
    bp[b]=pb;
    barea[b] += barea[b1]; barea[b1]=0.0;
    bforgetindex(b1);
    i=ii=iseed;  /* ii marks the endpoint */
    do { cadj[i][0]=b; i=vnbr[i][1];
    } while (i!=ii);
    bordpop(iseed,pressuredelta/2.0);
  }
  else {  /* case of 'b==b1' */
    /* for now, this case does *not* arise */
  }
  /*
   * Record edge losses for information purposes...
   */
  elosscount++;
  return;  /* Successful exit */
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  B P I N C H ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	bpinch(short i)
 *
 *	Arguments:	i	= index of vertex lying on the Plateau border
 *
 *	Return value:	none
 *
 *	Action:		This subroutine enacts Plateau border pinching, which
 *			is one of the two basic topological changes.  The
 *			first step in the algorithm is to fill the
 *			`wheel[][]' array which is essentially a list of the
 *			points lying around the perimeter of the Plateau
 *			border.  Then we loop around all possible pairs of
 *			arcs on the perimeter of the border which are at
 *			least *two* steps away from each other (therefore it
 *			follows that a 3-sided border never pinches, since
 *			each arc is adjacent to the other two).  Some
 *			elementary geometry is used to check if these two arcs
 *			overlap or touch each other.  If so, then a `pinch'
 *			is detected and the Plateau border gets split into
 *			two parts.  Finally, `bpinch()' is called on each of
 *			the two new parts in case further `pinching' is
 *			possible.
 *
 *****************************************************************************/
void bpinch(i)
short i;
{
  /* Make a wheel around the perimeter of the border. The period indices
  are all changed to be relative to 'i'=wheel[0][0]. */
  short b, b1, bb, bi, bj, i1, j, j1, ii, ii1, jj, jj1, k, kk, kk1, c1, c2, je,
        l, m, maxl;
  short nwh, wheel[MWHEEL][2], perconcat(), vgetindex(), bgetindex(),
        bsides();
  REAL x1, y1, x11, y11, x2, y2, x21, y21, alpha1, alpha2, p1, p2, pb,
       r1, r2, r, f, d, dx, dy, xc1, yc1, xc2, yc2, th, th1, ba, ba1, oba;
  REAL linlen(), barcangle(), linangle(), pimod(), twopimod(), bordarea(), quadarea();
  boolean la11, la12, la21, la22, straight1, straight2, flag, arcbrk;
  boolean arccentre(), vincrement(), larc();
  void bpinch(), setlarc(), trans(), putinbox(), bordpop();
  /*
   * Initialise the objects we are operating upon and, in particular,
   * we create a `wheel[][]' of points containing all the points around
   * the perimeter of the Plateau border.  The matrix wheel lists the
   * points in the format:
   * (vertex index = wheel[i][0], periodic index = wheel[i][1] )
   * and the points are listed in *anti-clockwise* order as we move
   * around the border perimeter.
   */
  b=cadj[i][0]; pb=bp[b];
  nwh=0;
  wheel[0][0]=i;
  i1=vper[i][1]; wheel[0][1]=i1 & LARC;
  nwh++;
  je=i;
  while ((j=vnbr[i][1])!=je) {
    j1=perconcat(i1,vper[j][1]); setlarc(&j1,larc(j,1));
    wheel[nwh][0]=j; wheel[nwh][1]=i1; nwh++;
    i=j; i1=j1;
  }
  /*
   * Now loop around every arc on the wheel to test for pinching...
   */
  maxl=nwh;
  for (l=0; l<maxl; l++) {
    i=wheel[l][0]; i1=wheel[l][1];
    ii=wheel[(l+1)%nwh][0]; ii1=wheel[(l+1)%nwh][1];
    trans(vx[i],vy[i],i1,&x1,&y1);
    trans(vx[ii],vy[ii],ii1,&x11,&y11);
    la11=larc(i,1);
    alpha1=barcangle(x1,y1,x11,y11,p1=cp[cadj[i][2]],pb,la11,&arcbrk);
    straight1= !arccentre(x1,y1,x11,y11,alpha1,&xc1,&yc1);
    /*
     * Loop around all other arcs which are at least two steps further
     * away along the perimeter (which guarantees that the pairs of
     * arcs are `opposing' arcs).
     */
    for (m=(l+2)%nwh; (l-m+nwh)%nwh>=2; m= (m+1)%nwh) {
      j=wheel[m][0]; j1=wheel[m][1];
      jj=wheel[(m+1)%nwh][0]; jj1=wheel[(m+1)%nwh][1];
      trans(vx[j],vy[j],j1,&x2,&y2);
      trans(vx[jj],vy[jj],jj1,&x21,&y21);
      la21=larc(j,1);
      alpha2=barcangle(x2,y2,x21,y21,p2=cp[cadj[j][2]],pb,la21,&arcbrk);
      straight2= !arccentre(x2,y2,x21,y21,alpha2,&xc2,&yc2);
      flag=FALSE;
      th=linangle(xc1,yc1,xc2,yc2); th1=linangle(xc1,yc1,x11,y11);
      /*
       * Test whether the angle of the first arc intersects the line
       * joining the centres of the two arcs (this being a necessary
       * but not a sufficient condition for the pinching of the 
       * two arcs).
       * The form of this test relies on the fact that the order of
       * points around the wheel is anti-clockwise.
       */
      if (twopimod(th-th1)<alpha1) {
        la12= (twopimod(alpha1+th1-th)>PI);
        la21= (twopimod(th-th1)>PI);
        th=twopimod(th+PI);
        th1=linangle(xc2,yc2,x21,y21);
        /*
         * Test whether the angle of the *second* arc intersects the line
         * joining the centres of the two arcs.
         * The form of this test relies on the fact that the order of
         * points around the wheel is anti-clockwise.
         */
        if (twopimod(th-th1)<alpha2) {
          la22= (twopimod(alpha2+th1-th)>PI);
          la11= (twopimod(th-th1)>PI);
	  /*
	   * Rule out the (perhaps remote) possibility that the
	   * Plateau border is twisted around so that both arcs
	   * are pointing in the same direction.
	   */
          if (quadarea(x1,y1,x11,y11,x2,y2,x21,y21) > 0.0) {
	    /*
	     * Having already checked that both arcs intersect the
	     * line joining the centres of the arcs, it is sufficient
	     * to show that the sum of the arc radii is larger than
	     * the distance between arc centres.
	     */
            r=linlen(xc1,yc1,xc2,yc2);
            if (!straight1) {
              if (!straight2) {
                /* General case */
                r1=BRADIUS(p1,pb); r2=BRADIUS(p2,pb);
                if ((r1+r2)+filmwid>r) {
                  flag=TRUE;
		  /*
		   * `f' tells you where the pinching arcs touch
		   * or, more realistically, half way between where
		   * the two arcs overlap -- expressed as a fraction
		   * of the vector (xc2-xc1, yc2-yc1).
		   */
                  f=r1+(r-(r1+r2))/2.0; f /= r;
                }
              }
              else {
                /* In this case, 'arccentre()' is assumed to have returned
                 * the midpoint of the line segment */
                r1=BRADIUS(p1,pb);
                if (r1+filmwid>r) {
                  flag=TRUE;
                  f=r1+(r-r1)/2.0; f /= r;
                }
              }
            }
            else {
              if (!straight2) {
                r2=BRADIUS(p2,pb);
                if (r2+filmwid>r) {
                  flag=TRUE;
                  f=(r-r2)/2.0; f /= r;
                }
              }
              else {
                if (filmwid>r) { flag=TRUE; f=0.5; }
              }
            }
          }
          if (flag) {   /* pinch found */
            bpinchcount++;
            /*
	     * First enact the topological change and
	     * fix the endpoints of the new edge.
	     */
            k=vgetindex();
            kk=vgetindex();
            b1=bgetindex();
	    /* (xc1+dx, yc1+dx) is where the arcs touch */
            dx=f*(xc2-xc1); dy=f*(yc2-yc1);
            d=sqrt(dx*dx+dy*dy);
            vx[k]=xc1+dx-5.0*minvvlen*dy/d;
            vy[k]=yc1+dy+5.0*minvvlen*dx/d;
            vx[kk]=xc1+dx+5.0*minvvlen*dy/d;
            vy[kk]=yc1+dy-5.0*minvvlen*dx/d;
            dvx[k]=dvy[k]=dvx[kk]=dvy[kk]=0.0;
            vnbr[k][0]=kk; vnbr[k][1]=jj; vnbr[k][2]=i;
            vnbr[kk][0]=k; vnbr[kk][1]=ii; vnbr[kk][2]=j;
            vnbr[i][1]=k; vnbr[jj][2]=k;
            vnbr[j][1]=kk; vnbr[ii][2]=kk;
            c1=cadj[i][2]; c2=cadj[j][2];
            cadj[k][0]=b; cadj[k][1]=c1; cadj[k][2]=c2;
            cadj[kk][0]=b1; cadj[kk][1]=c2; cadj[kk][2]=c1;
            /*
	     * Since `b1' is a newly created Plateau border, it
	     * is necessary to tell all the points on its perimeter
	     * about its existence.
	     */
            bj=je=kk;
            do { cadj[bj][0]=b1; bj=vnbr[bj][1];
            } while (bj!=je);
            /*
	     * Now update the array 'vper', noting that the periodic
             * indices i1,j1,ii1,jj1 were defined so that they translate
	     * their points into the vicinity of the first point on the
	     * wheel (wheel[0][0]).
	     */
            vper[k][0]=0; vper[k][1]=jj1 & PERMASK; vper[k][2]=i1 & PERMASK;
            setlarc(&vper[k][1],la11); setlarc(&vper[k][2],la12);
            vper[kk][0]=0; vper[kk][1]=ii1 & PERMASK; vper[kk][2]=j1 & PERMASK;
            setlarc(&vper[kk][1],la21); setlarc(&vper[kk][2],la22);
            vper[i][1]=PERFN(-PERX(i1),-PERY(i1));
            setlarc(&vper[i][1],la12);
            vper[jj][2]=PERFN(-PERX(jj1),-PERY(jj1));
            setlarc(&vper[jj][2],la11);
            vper[j][1]=PERFN(-PERX(j1),-PERY(j1));
            setlarc(&vper[j][1],la22);
            vper[ii][2]=PERFN(-PERX(ii1),-PERY(ii1));
            setlarc(&vper[ii][2],la21);
            /* This clinches the correct values of vper for 'k', 'kk'
             * and their nbring vertices */
            putinbox(k); putinbox(kk);
	    /*
	     * Now update the remaining network information...
	     */
            nbsides[b]=bsides(b); nbsides[b1]=bsides(b1);
            (ncsides[c1])++; (ncsides[c2])++;
            bp[b1]=pb;
            bordpop(k,pressuredelta/2.0); bordpop(kk,pressuredelta/2.0);
            oba=barea[b];
            ba=bordarea(k,&arcbrk); ba1=bordarea(kk,&arcbrk);
            barea[b]=oba*(ba/(ba+ba1)); barea[b1]=oba*(ba1/(ba+ba1));
	    /*
	     * Call `bpinch()' recursively on each of the newly created
	     * Plateau borders.
	     */
            if (nbsides[b]>3) bpinch(k);
            if (b1<b && nbsides[b1]>3) bpinch(kk);
            return;
          }
        }
      }
    }
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *    V I N C R E M E N T ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	boolean vincrement(short i, REAL dx, dy)
 *
 *	Arguments:	i	= index of vertex to increment
 *			(dx, dy) = size of vertex increment to attempt
 *
 *	Return value:	none
 *
 *	Action:		This subroutine achieves what you might think should
 *			be a trivial task -- it adds (dx, dy) the vertex
 *			coordinates (vx[i], vy[i]).  However making this
 *			increment is slightly hazardous since it is possible
 *			to move the vertex into a position where the network
 *			becomes geometrically impossible(!).  Therefore this
 *			routine institutes a number of checks to make sure the
 *			geometry remains intact.  If a vertex is moved so far
 *			away that a border arc no longer fits between itself
 *			and its neighbour than things are patched up by
 *			shifting this vertex a bit.
 *
 *****************************************************************************/
boolean vincrement(i,dx,dy)
short i;
REAL dx,dy;
{
  short j, j1, jj, jj1, ii, k, kk, c, cc;
  REAL x1, y1, x11, y11, x12, y12, d, dd, r, rr, dx1, dy1, f, th, dxx, dyy,
       ta, tb, tc, pb;
  REAL linlen(); 
  void setlarc(), trans(), putinbox(), vnbrxy();
  boolean la, found, larc(), vincrement();
  if (dx==0.0 && dy==0.0) return TRUE;
  else {
    /*
     * In fact we are being a little illogical here, since we are
     * assuming that the increment passed as (dx, dy) is
     * equal to (dvx[i], dvy[i]) !
     *
     * It is important to subtract the increment (dx, dy) from
     * the increment arrays (dvx[i], dvy[i]) as soon as they are
     * added to (vx[i], vy[i]).  Otherwise, on account of this routine
     * calling itself, the vertex could get incremented more than
     * once.
     */
    vx[i] += dx; vy[i] += dy;
    dvx[i] -= dx; dvy[i] -= dy;
    /*
     * Look at each of the two neighbouring arcs to see if they are
     * overstretched when vertex 'i' is incremented...
     */
    x1=vx[i]; y1=vy[i]; pb=bp[cadj[i][0]];
    for (k=1; k<=2; k++) {
      kk=k%2+1;
      j=vnbr[i][k]; j1=vper[i][k];
      la=larc(i,k);
      trans(vx[j],vy[j],j1,&x11,&y11);
      d=linlen(x1,y1,x11,y11);
      c=cadj[i][kk]; r=BRADIUS(cp[c],pb);
      /*
       * The largest possible distance between two neighbouring
       * points is given by the largest chord of the arc joining
       * them  = 2.0*r.  In practice we don't let arcs extend to their
       * full diameter (it is numerically useful to have some play left)
       * so `cosminang', which is slightly less than one, is a parameter
       * which limits the largest size of arc.
       */
      if (d>(ta=2.0*r*cosminang)) {
	/*
	 * If the arc is ``broken'' then the first possible remedy is to
	 * give the neighbouring point the increment which is due to
	 * it in any case from (dvx[j], dvy[j]).
	 */
        if (vincrement(j,dvx[j],dvy[j])) {
          vnbrxy(i,k,&x11,&y11);
          d=linlen(x1,y1,x11,y11);
        }
        if (d>(ta=2.0*r*cosminang)) {
	  /*
	   * If that didn't work then we reach here.  The ultimate
	   * solution (not very elegant but it works) is to move
	   * vertex `i' just a little bit closer to its neighbour.
	   */
          f=1.0-2.0*r*cosminang/d; f *= 1.1;
          dx1=f*(x11-x1); dy1=f*(y11-y1);
          vx[i] += dx1; vy[i] += dy1;
          /*
	   * Now check whether the other border arc is now out of place
	   */
          jj=vnbr[i][kk]; jj1=vper[i][kk];
          trans(vx[jj],vy[jj],jj1,&x12,&y12);
          dd=linlen(x1,y1,x12,y12);
          cc=cadj[i][k]; rr=BRADIUS(cp[cc],pb);
          if (dd>(tb=2.0*rr*cosminang)) {
	    /*
	     * If the arc is ``broken'' then the first possible remedy is to
	     * give the neighbouring point the increment which is due to
	     * it in any case from (dvx[jj], dvy[jj]).
	     */
            if (vincrement(jj,dvx[jj],dvy[jj])) {
              vnbrxy(i,kk,&x12,&y12);
              dd=linlen(x1,y1,x12,y12);
            }
            if (dd>(tb=2.0*rr*cosminang)) {
	      /*
	       * If that didn't work then we reach here.  The solution
	       * here is a little bit trickier than last time because
	       * in fixing *this* arc, we don't want to end up
	       * breaking the previous one:
	       *
               * The point 'i' is moved to the apex of a triangle of
               * side lengths 'ta' and 'tb'
	       */
              tc=linlen(x11,y11,x12,y12);
              dxx=(x12-x11)/tc; dyy=(y12-y11)/tc;
              th=acos((ta*ta+tc*tc-tb*tb)/(2.0*ta*tc));
              vx[i]=vx[j]+ta*(cos(th)*dxx-sin(th)*dyy);
              vy[i]=vy[j]+ta*(sin(th)*dxx+cos(th)*dyy);
            }
            else if (nbsides[cadj[i][0]]>3) {
              /* if successfully de-stretched and not a 3-border then 'pop' */
              setlarc(&vper[i][k],!la);
              setlarc(&vper[j][k%2+1],!la);
            }
          }
        }
      }
    }
  }
  return TRUE;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  H E N C K Y ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	hencky(eps)
 *
 *	Arguments:	eps	= Hencky strain parameter
 *
 *	Return value:	none
 *
 *	Action:		Apply Hencky strain to the network.  This has the
 *			effect of multiplying the height of the box by
 *			`exp(eps)' and multiplying the width of the box by
 *			`exp(-eps)', thus preserving the area of the network.
 *			The vertices are all transformed proportionately so
 *			that they lie in roughly the right positions.
 *			Equilibration is *not* carried out as part of this
 *			step.
 *
 *****************************************************************************/
void hencky(eps)
REAL eps;
{
  short ii, i, j, k, c;
  boolean found[MCELL], vincrement(), arcbrk;
  REAL x1, y1, x2, y2, eofeps;
  REAL equil(), cellarea();
  void polyplot(), plerror(), foamplot();
  /*
   * Strain the periodic box.
   */
  henckyeps += eps;
  eofeps=exp(eps);
  boxwid /= eofeps; boxhgt *= eofeps;
  /*
   * Perform a simple linear transformation of the
   * vertices to put them in approximately the
   * right positions...
   */
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    x2=(x1=vx[i])/eofeps; y2=(y1=vy[i])*eofeps;
    dvx[i]=x2-x1; dvy[i]=y2-y1;
  }
  /*
   * Call equil with an `adjust' flag, so that *no*
   * actual equilibration is done, just some checks
   * when the increments (dvx[], dvy[]) are added.
   */
  equil(TRUE);
}
