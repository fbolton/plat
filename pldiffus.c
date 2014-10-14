#include "include.h"

/******************************/
/* Diffusion related routines */
/******************************/

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  D I F F U S E ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	diffuse()
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		Implements intercellular diffusion of gas using a
 *			modified von Neumann law.  Essentially it is assumed
 *			that gas can diffuse through that part of a cell-cell
 *			arc which is not covered by a Plateau border.
 *
 *			The rate of diffusion is determined by `diffuserate'
 *			and this is scaled to give a rate of diffusion
 *			proportional to `Abar' (the average area per cell).
 *			This means that cell areas change by a certain
 *			percentage per step.
 *
 *			NB:  This routine is very simple but it doesn't do
 *			everything it needs to do.  In particular when cells
 *			get very small as a consequence of diffusion, problems
 *			occur (which are not solved in the current program).
 *			This subroutine has not been used much...
 *
 *****************************************************************************/
void diffuse()
{
  short i, ii, j, c1, c2;
  REAL abar, da, dt, x1, y1, x2, y2, p1, p2;
  REAL carclen();
  void vnbrxy();
  abar=boxwid*boxhgt/((REAL) nc);
  dt=diffuserate*abar;
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    j=vnbr[i][0];
    if (j>i) {
      x1=vx[i]; y1=vy[i];
      vnbrxy(i,0,&x2,&y2);
      p1=cp[c1=cadj[i][1]]; p2=cp[c2=cadj[i][2]];
      da=dt*(p1-p2)*carclen(x1,y1,x2,y2,p1,p2);
      darea[c2] += da; darea[c1] -= da;
    }
  }
  tfoam += dt;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *  M A K E F R A C T I O N ( )     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	makefraction(REAL phi0, boolean plotflag)
 *
 *	Arguments:	phi0	= target value for the new area fraction `phi'
 *			plotflag = TRUE selects plotting at each iteration
 *
 *	Return value:	none
 *
 *	Action:		This routine changes the area fraction of a wet foam
 *			until it reaches the target fraction `phi0' (which
 *			may be either greater or lesser than the existing 
 *			fraction `phi').  Since the initialisation routines
 *			can only generate foams with very small Plateau
 *			borders (i.e. `phi' close to 0.99) this subroutine
 *			is crucially important for generating wetter foams.
 *
 *			However there is currently a serious bug in this
 *			algorithm which means that it does not cope well with
 *			very wet foams (this ought to be fixed soon).
 *
 *			The outline of the algorithm is straightforward and
 *			goes as follows:
 *			(i) Firstly choose a step size `dphi' by which the
 *			current area fraction will be changed (we do not want
 *			to change `phi' all in one go).
 *			(ii) Then change the target areas of all the cells so
 *			as to match the new area fraction.
 *			(iii) Call `phirelax(phi)' so as to change the
 *			geometry of the network to match the new area fraction.
 *			The routine `phirelax()' works by changing the
 *			curvature of the border arcs (thereby changing the
 *			areas of adjacent cells by a small amount).
 *			(iv) Then equilibrate the vertices of the network by
 *			calling `equil()' (which works by imposing the
 *			condition that arcs match up smoothly).
 *			(v) Iterate until the target area fraction `phi0'
 *			is reached.
 *
 *****************************************************************************/
void makefraction(phi0, plotflag)
REAL phi0;
boolean plotflag;
{
  short i, ii, j;
  REAL phi, dphi, atot;
  REAL equil();
  void foamplot(), setconstants(), phirelax();
#ifdef IRIS
  void zoom();
#endif
  /*
   * Calculate `atot' which is what the total area of gas *ought* to be
   * at present (even though geometrically the network may have cells which
   * are slightly the wrong size).
   *
   * Note that `darea[i]' records the difference between the area which a
   * cell `i' *ought* to have, and the area which (geometrically) it actually
   * contains.
   */
  atot=0.0;
  for (ii=0; ii<nc; ii++) {
    i=clist[ii];
    atot += carea[i]+darea[i];
  }
  for (ii=0; ii<nbub; ii++) {
    i=bublist[ii];
    atot += carea[i]+darea[i];
  }
  phi=(atot)/(boxwid*boxhgt);
  /*
   * Down to business...
   */
  if (phi>phi0) {
    /*
     * Case 1:  Shrink the area fraction down to `phi0'...
     */
    /*
     * Aim to change the area fraction in steps `dphi' which
     * must not be overly large (because `phirelax()' cannot
     * cope with large steps).
     */
    dphi=(1.0-phi)*PHIFRACTION;
    if (fabs(phi-phi0)<fabs(dphi)) dphi=phi-phi0;
    while (fabs(dphi/phi) > areasup) {
      /* Note: areasup = `negligibly small area' */
      /*
       * Change the target areas of all the cells and recalculate
       * `atot' and `phi' ( = current target for area fraction)
       */
      atot=0.0;
      for (ii=0; ii<nc; ii++) {
        i=clist[ii];
        darea[i] -= (carea[i]+darea[i])*dphi/phi;
        atot += carea[i]+darea[i];
      }
      for (ii=0; ii<nbub; ii++) {
        i=bublist[ii];
        darea[i] -= (carea[i]+darea[i])*dphi/phi;
        atot += carea[i]+darea[i];
      }
      phi=(atot)/(boxwid*boxhgt);
      /*
       * Actually change the area fraction of the network
       * (by tweaking the curvature of the Plateau border arcs)
       *
       * ...and do it to a pretty high degree of accuracy (i.e. the parameter
       * `areasup' is very small [or certainly ought to be] ) otherwise
       * the equilibration would bomb out.
       */
      phirelax(phi);
      /*
       * Call `equil()' to relax the position of the vertices...
       */
      j=0;
      while (equil(0)>equilsup && ((REAL) j)<smaxiter) {
        setconstants();
        if (plotflag) {
#ifdef IRIS
          zoom();
#endif
          foamplot(0);
        }
        j++;
      }
      /*
       * Get a new step size `dphi' before continuing the loop...
       */
      dphi=(1.0-phi)*PHIFRACTION;
      if (fabs(phi-phi0)<fabs(dphi)) dphi=phi-phi0;
    }
  }
  else {
    /*
     * Case 2:  Increase the area fraction up to `phi0'...
     */
    /*
     * Aim to change the area fraction in steps `dphi' which
     * must not be overly large (because `phirelax()' cannot
     * cope with large steps).
     */
    dphi=(1.0-phi)*PHIFRACTION;
    if (fabs(phi-phi0)<fabs(dphi)) dphi=phi0-phi;
    while (fabs(dphi/phi) > areasup) {
      /*
       * Change the target areas of all the cells and recalculate
       * `atot' and `phi' ( = current target for area fraction)
       */
      atot=0.0;
      for (ii=0; ii<nc; ii++) {
        i=clist[ii];
        darea[i] += (carea[i]+darea[i])*dphi/phi;
        atot += carea[i]+darea[i];
      }
      for (ii=0; ii<nbub; ii++) {
        i=bublist[ii];
        darea[i] += (carea[i]+darea[i])*dphi/phi;
        atot += carea[i]+darea[i];
      }
      phi=(atot)/(boxwid*boxhgt);
      /*
       * Actually change the area fraction of the network
       * (by tweaking the curvature of the Plateau border arcs)
       *
       * ...and do it to a pretty high degree of accuracy (i.e. the parameter
       * `areasup' is very small [or certainly ought to be] ) otherwise
       * the equilibration would bomb out.
       */
      phirelax(phi);
      /*
       * Call `equil()' to relax the position of the vertices...
       */
      j=0;
      while (equil(0)>equilsup && ((REAL) j)<smaxiter) {
        setconstants();
        if (plotflag) {
#ifdef IRIS
          zoom();
#endif
          foamplot(0);
        }
        j++;
      }
      /*
       * Get a new step size `dphi' before continuing the loop...
       */
      dphi=(1.0-phi)*PHIFRACTION;
      if (fabs(phi-phi0)<fabs(dphi)) dphi=phi0-phi;
    }
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *P H I R E L A X ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	phirelax(REAL phi)
 *
 *	Arguments:	phi	= target value of the area fraction
 *
 *	Return value:	none
 *
 *	Action:		Relax the network till it reaches a volume fraction
 *			given by the target value `phi' (this is a bit of a
 *			misnomer since the network is not terribly relaxed by
 *			the end of this routine).
 *
 *			It operates by changing the pressure of the Plateau
 *			borders, thereby effectively changing the area in the
 *			cells.  Note that since the cell areas are changed
 *			merely by tweaking the curvature of the Plateau border
 *			arcs, this routine is only able to change `phi' by a
 *			fairly small amount.
 *
 *****************************************************************************/
void phirelax(phi)
REAL phi;
{
  short i, ii, b;
  REAL phi1, phi2, dphi, dphidbp, delp, delp2, dbp, r;
  boolean found[MBORD], arcbrk;
  void incbp(), bordpop();
  REAL phifn(), fsign(), bordarea();
  delp=pressuredelta; delp2=delp/2.0;
  dphi=phi-phifn();
  /*
   * Essentially this loop tries to alter the area fraction
   * so that it equals `phi'.  It does this by calculating
   * the derivative \delta\phi / \delta bp and using this to
   * guess an amount by which the border pressure (= bp) must
   * be changed.
   *
   * It's a wee bit crude (in fact for *large* Plateau borders it is not all
   * that successful either --> significant bug (I promise I'll fix it!)  )
   */
  while (fabs(dphi/phi) > areasup) {
    incbp(-delp2); phi1=phifn();
    incbp(delp); phi2=phifn();
    incbp(-delp2);
    /* Possibly include damping for this step =0.5? */
    dphidbp=0.5*(phi2-phi1)/delp;
    dbp= PHIDAMP*dphi/dphidbp;
    if (fabs(dbp)>fabs(0.1*bpav)) dbp=fsign(dbp)*fabs(0.1*bpav);
    incbp(dbp);
    dphi=phi-phifn();
  }
  /* correct the bubble areas */
  for (ii=0; ii<nbub; ii++) {
    i=bublist[ii];
    r=BRADIUS(cp[i],bpav);
    carea[i]= PI*r*r; darea[i]=0.0;
  }
  /* correct the border areas */
  for (i=0; i<onb; i++) found[i]=FALSE;
  for (ii=0; ii<nb; ii++) {
    i=vlist[ii];
    if (!found[b=cadj[i][0]]) {
      found[b]=TRUE;
      barea[b]=bordarea(i,&arcbrk);
      if (arcbrk) {
        bordpop(i,pressuredelta/2.0);
        barea[b]=bordarea(i,&arcbrk);
      }
    }
  }
  /* ...but don't worry about correcting the cell areas,
   * that business is taken care of in the routine `makefraction()'
   */
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *    I N C B P ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	incbp(REAL dbp)
 *
 *	Arguments:	dbp	= increment to the Plateau border pressure
 *
 *	Return value:	none
 *
 *	Action:		Increments the Plateau border pressure by an amount
 *			`dbp'.
 *
 *****************************************************************************/
void incbp(dbp)
REAL dbp;
{
  short i, ii;
  for (ii=0; ii<nb; ii++) {
    i=blist[ii];
    bp[i] += dbp;
  }
  bpav += dbp;
  /*
   * This is a tricky point:
   *
   * In order to avoid implementing isolated bubbles, the program allows
   * single sided cells to exist (i.e. the cells are almost a complete circle)
   * which just `kiss' their neighbours.  These cells are *not* able to change
   * their curvature in reaction to a change in Plateau border pressure so
   * their internal pressure has to track changes in Plateau border pressure.
   */
  for (ii=0; ii<nc; ii++) {
    i=clist[ii];
    if (ncsides[i]==1) cp[i] += dbp;
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *    P H I F N ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL phifn()
 *
 *	Arguments:	none
 *
 *	Return value:	Returns the area fraction \phi as calculated explicitly
 *
 *	Action:		(same as return value)
 *
 *****************************************************************************/
REAL phifn()
{
  short i, ii, k, c;
  REAL atot, r, cellarea(), acell;
  boolean found[MCELL], arcbrk;
  void cellpop();
  /*
   * Calculate the area fraction by adding together the areas of all of the
   * cells (calculated explicitly) and then dividing by the total area of the
   * network (which is given by boxwid*boxhgt).
   */
  atot=0.0;
  for (ii=0; ii<nc; ii++) { found[ii]=FALSE; }
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    for (k=1; k<3; k++) {
      if (!found[c=cadj[i][k]]) {
        found[c]=TRUE;
        acell=cellarea(i,k,&arcbrk);
        if (arcbrk) {
          cellpop(i,k,pressuredelta);
          acell=cellarea(i,k,&arcbrk);
        }
        atot += acell;
      }
    }
  }
  for (ii=0; ii<nbub; ii++) {
    i=bublist[ii];
    r=BRADIUS(cp[i],bpav);
    atot += PI*r*r;
  }
  return (atot/(boxwid*boxhgt));
}
