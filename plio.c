#include "include.h"

/****************/
/* I/O routines */
/****************/

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *V F O A M O U T ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	vfoamout(char *s)
 *
 *	Arguments:	s[]	= output filename
 *
 *	Return value:	none
 *
 *	Action:		Writes the Voronoi network the file whose name is
 *			given by the argument s[].
 *
 *****************************************************************************/
void vfoamout(s)
char *s;
/* The preferred name is "vfoam.pb" */
{
  FILE *fp;
  short i, j;
  if (! (fp=fopen(s,"w")) ) {
    plerror("error: vfoamout: couldn't open file to write");
    return;
  }
  fprintf(fp,"%d %d\n",nv,nc);
  for (i=0; i<nv; i++) {
    fprintf(fp,"%f %f\n",vx[i],vy[i]);
    for (j=0; j<3; j++) fprintf(fp,"%d ",vorvnbr[i][j]);
    fprintf(fp,"\n");  
    for (j=0; j<3; j++) fprintf(fp,"%d ",vorcadj[i][j]);
    fprintf(fp,"\n");  
    for (j=0; j<3; j++) fprintf(fp,"%d ",vorvper[i][j]);
    fprintf(fp,"\n");  
  }
  fclose(fp);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     * V F O A M I N ( )    *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	vfoamin(char *s)
 *
 *	Arguments:	s[]	= name of file to read
 *
 *	Return value:	none
 *
 *	Action:		Reads in details of a Voronoi network from the file
 *			named in argument s[].
 *
 *****************************************************************************/
void vfoamin(s)
char *s;
/* The preferred name is "vfoam.pb" */
{
  FILE *fp;
  short i, j;
  if (! (fp=fopen(s,"r")) ) {
    plerror("error: vfoamin: couldn't open file to read");
    return;
  }
  fscanf(fp,"%d %d",&nv,&nc);
  onv=nv; onc=nc;
  for (i=0; i<nv; i++) {
    fscanf(fp,"%f %f",&vx[i],&vy[i]);
    for (j=0; j<3; j++) fscanf(fp,"%d",&vorvnbr[i][j]);
    for (j=0; j<3; j++) fscanf(fp,"%d",&vorcadj[i][j]);
    for (j=0; j<3; j++) fscanf(fp,"%d",&vorvper[i][j]);
  }
  fclose(fp);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     * F O A M O U T ( )    *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	foamout(char *s)
 *
 *	Arguments:	s[]	= name of output file
 *
 *	Return value:	none
 *
 *	Action:		Writes a Plateau topology network to the file named
 *			in the argument s[].
 *
 *****************************************************************************/
void foamout(char *s)
/* The preferred name is "foam.pb" */
{
  short i, ii, j;
  FILE *fp;
  if (! (fp=fopen(s,"w")) ) {
    plerror("error: foamout: couldn't open file to write");
    return;
  }
  fprintf(fp,"%d %d %d %d\n",nv,nc,nb,nbub);
  fprintf(fp,"%d %d %d\n",onv,onc,onb);
  fprintf(fp,"%f %f %f %f\n",boxwid,boxhgt,bpav,henckyeps);
  for (ii=0; ii<nv; ii++) {
    fprintf(fp,"%d\n",i=vlist[ii]);
    fprintf(fp,"%f %f\n",vx[i],vy[i]);
    for (j=0; j<3; j++) fprintf(fp,"%d ",vnbr[i][j]);
    fprintf(fp,"\n");
    for (j=0; j<3; j++) fprintf(fp,"%d ",vper[i][j]);
    fprintf(fp,"\n");
    for (j=0; j<3; j++) fprintf(fp,"%d ",cadj[i][j]);
    fprintf(fp,"\n");
  }
  for (ii=0; ii<nc; ii++) {
    fprintf(fp,"%d\n",i=clist[ii]);
    fprintf(fp,"%d\n",ncsides[i]);
    fprintf(fp,"%f\n",cp[i]);
    fprintf(fp,"%f\n",carea[i]+darea[i]);
  }
  for (ii=0; ii<nb; ii++) {
    fprintf(fp,"%d\n",i=blist[ii]);
    fprintf(fp,"%d\n",nbsides[i]);
    fprintf(fp,"%f\n",bp[i]);
    fprintf(fp,"%f\n",barea[i]);
  }
  for (ii=0; ii<nbub; ii++) {
    fprintf(fp,"%d\n",i=bublist[ii]);
    fprintf(fp,"%f\n",cp[i]);
    fprintf(fp,"%f\n",carea[i]+darea[i]);
    fprintf(fp,"%f %f\n",cx[i],cy[i]);
  }
  fclose(fp);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  F O A M I N ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	foamin(char *s)
 *
 *	Arguments:	s[]	= name of input file
 *
 *	Return value:	none
 *
 *	Action:		Inputs a Plateau topology network from the file name
 *			given by argument s[].
 *
 *****************************************************************************/
void foamin(char *s)
/* The preferred name is "foam.pb" */
{
  FILE *fp;
  short i, ii, j;
  int d1, d2, d3, d4;
  float f1, f2, f3, f4;
  if (! (fp=fopen(s,"r")) ) {
    plerror("error: foamin: couldn't open file to read");
    return;
  }
  fscanf(fp,"%d %d %d %d",&d1,&d2,&d3,&d4);
  nv=(short) d1; nc=(short) d2; nb=(short) d3; nbub=(short) d4;
  fscanf(fp,"%d %d %d",&d1,&d2,&d3);
  onv=(short) d1; onc=(short) d2; onb=(short) d3;
  fscanf(fp,"%f %f %f %f",&f1,&f2,&f3,&f4);
  boxwid=(REAL) f1; boxhgt=(REAL) f2; bpav=(REAL) f3; henckyeps=(REAL) f4;
  for (ii=0; ii<nv; ii++) {
    fscanf(fp,"%d",&d1); i=(short) d1; vlist[ii]=i;
    fscanf(fp,"%f %f",&f1,&f2);
    vx[i]=(REAL) f1; vy[i]=(REAL) f2;
    for (j=0; j<3; j++) { fscanf(fp,"%d",&d1); vnbr[i][j]=(short) d1; }
    for (j=0; j<3; j++) { fscanf(fp,"%d",&d1); vper[i][j]=(short) d1; }
    for (j=0; j<3; j++) { fscanf(fp,"%d",&d1); cadj[i][j]=(short) d1; }
  }
  for (ii=0; ii<nc; ii++) {
    fscanf(fp,"%d",&d1); i=(short) d1; clist[ii]=i;
    fscanf(fp,"%d",&d1); ncsides[i]=(short) d1;
    fscanf(fp,"%f",&f1); cp[i]=(REAL) f1;
    fscanf(fp,"%f",&f1); carea[i]=(REAL) f1;
    darea[i]=0.0;
  }
  for (ii=0; ii<nb; ii++) {
    fscanf(fp,"%d",&d1); i=(short) d1; blist[ii]=i;
    fscanf(fp,"%d",&d1); nbsides[i]=(short) d1;
    fscanf(fp,"%f",&f1); bp[i]=(REAL) f1;
    fscanf(fp,"%f",&f1); barea[i]=(REAL) f1;
  }
  for (ii=0; ii<nbub; ii++) {
    fscanf(fp,"%d",&d1); i=(short) d1; bublist[ii]=i;
    fscanf(fp,"%f",&f1); cp[i]=(REAL) f1;
    fscanf(fp,"%f",&f1); carea[i]=(REAL) f1;
    fscanf(fp,"%f %f",&f1,&f2); cx[i]=(REAL) f1; cy[i]=(REAL) f2;
    darea[i]=0.0;
  }
  fclose(fp);
}

#ifndef DELAUNEY
/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     V F O A M P L O T ( )  *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	vfoamplot()   (Version I)
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		Plots a Voronoi network on the current graphic devices.
 *
 *****************************************************************************/
void vfoamplot()
{
  short i, ii, k, j, j1;
  REAL x1, y1, x2, y2;
#ifdef DEBUG
  REAL rmin;
#endif
  void mglinit(), mglline(), mglarc(), mglrect(), mglclose(), trans();
  mglinit();
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    x1=vx[i]; y1=vy[i];
    for (k=0; k<3; k++) {
      j=vorvnbr[i][k]; j1=vorvper[i][k];
      if ((j>i) || ((j1 & PERMASK)!=0)) {
        trans(vx[j],vy[j],j1,&x2,&y2);
        mglline(x1,y1,x2,y2);
      }
    }
  }
#ifdef DEBUG
  rmin=sqrt(boxwid*boxwid+boxhgt*boxhgt)*HARDRFRAC/2.0;
#endif
  k=0;
  for (ii=0; ii<nc; ii++) {
    i=clist[ii];
    x1=cx[i]; y1=cy[i];
    ptbox(&x1,&y1,&k,&k,&k);
    mglline(x1,y1,x1,y1);
#ifdef NOTNOW
    mglarc(x1,y1,rmin,0.0,2.0*PI);
#endif
  }
  mglrect(-0.5*boxwid,-0.5*boxhgt,0.5*boxwid,0.5*boxhgt);
  mglclose();
}
#else

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *      V F O A M P L O T ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	vfoamplot()   (Version II)
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		Plots the Delauney tiling of triangles to the current
 *			graphics devices.
 *
 *****************************************************************************/
void vfoamplot()
{
  short i, ii, k, j, j1, ca, cb, cc, ca1, cb1, cc1;
  REAL x1, y1, x2, y2, x3, y3;
#ifdef DEBUG
  REAL rmin;
#endif
  void mglinit(), mglline(), mglarc(), mglrect(), mglclose(), trans();
  mglinit();
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    ca=triang[i][0][0]; ca1=triang[i][0][1];
    cb=triang[i][1][0]; cb1=triang[i][1][1];
    cc=triang[i][2][0]; cc1=triang[i][2][1];
    trans(cx[ca],cy[ca],ca1,&x1,&y1);
    trans(cx[cb],cy[cb],cb1,&x2,&y2);
    trans(cx[cc],cy[cc],cc1,&x3,&y3);
    mglline(x1,y1,x2,y2); mglline(x2,y2,x3,y3); mglline(x3,y3,x1,y1);
  }
  mglrect(-0.5*boxwid,-0.5*boxhgt,0.5*boxwid,0.5*boxhgt);
  mglclose();
}
#endif

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *F O A M P L O T ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	foamplot(boolean hatchflag)
 *
 *	Arguments:	hatchflag = TRUE to switch on hatching of borders
 *
 *	Return value:	none
 *
 *	Action:		Plots the Plateau topology foam to the current
 *			graphics devices.  This routine is effectively a
 *			wrapper for the substantial routine `mglfoam()'.
 *
 *****************************************************************************/
void foamplot(boolean hatchflag)
{
  void mglinit(), mglfoam(), mglrect(), mglclose();
  mglinit();
  mglfoam(0.0,0.0,0.0,TRUE,TRUE,30.0,hatchflag);
  mglrect(-0.5*boxwid,-0.5*boxhgt,0.5*boxwid,0.5*boxhgt);
  mglclose();
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *P O L Y P L O T ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	polyplot(short mult, boolean hatchflag)
 *
 *	Arguments:	mult	= how many multiples of periodic box
 *			hatchflag = TRUE to switch on hatching of borders
 *
 *	Return value:	none
 *
 *	Action:		Like `foamplot()', except many periodic boxes are
 *			plotted.  The finished plot consists of `mult x mult'
 *			periodic boxes plotted to the current graphics devices.
 *
 *****************************************************************************/
void polyplot(short mult, boolean hatchflag)
{
  void mglinit(), mglfoam(), mglrect(), mglclose();
  REAL s, x, y, del, lowedge;
  del=1.0/((REAL) mult);
  s=del;
  lowedge= -0.5*(1.0-del);
  mglinit();
  for (x=lowedge; x<0.5; x += del)
    for (y=lowedge; y<0.5; y += del)
      mglfoam(x,y,s,(x==lowedge),(y==lowedge),30.0,hatchflag);
  mglrect(-0.5*boxwid,-0.5*boxhgt,0.5*boxwid,0.5*boxhgt);
  mglclose();
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *   C L U S T E R P L O T ( )*     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	clusterplot(short i)
 *
 *	Arguments:	i	= index of vertex lying on a Plateau border
 *
 *	Return value:	none
 *
 *	Action:		Plots the Plateau border, identified by the vertex `i'
 *			which it contains, to the current graphics devices.
 *			This subroutine is essentially a wrapper for the
 *			routine `mglbord()' which does most of the work.
 *
 *****************************************************************************/
void clusterplot(i)
short i;
{
  void vnbrxy(), mglinit(), mglline(), mglbord(), mglrect(), mglclose();
  REAL x1, y1, x2, y2;
  mglinit();
  x1=vx[i]; y1=vy[i]; vnbrxy(i,0,&x2,&y2);
  mglline(x1,y1,x2,y2);
  mglbord(i);
  mglrect(-0.5*boxwid,-0.5*boxhgt,0.5*boxwid,0.5*boxhgt);
  mglclose();
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *I N F O O P E N ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	infoopen()
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		Opens the file "info.pb" which will be used to record
 *			information and statistics about the foam in a
 *			tab delimited columns format.
 *
 *****************************************************************************/
static FILE *ff;

void infoopen()
{
  short i;
  void infotitles();
  if (! (ff=fopen("info.pb","w")) ) {
    plerror("error: infoopen: couldn't open file 'info.pb'");
    return;
  }
  fprintf(ff,"*\n");
  infotitles();
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *    I N F O T I T L E S ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	infotitles()
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		This routine prints the headings of the columns of
 *			information to the file "info.pb" (which must already
 *			have been opened by a call to `infoopen()' ).  The
 *			array `info_list[i]' contains TRUE or FALSE depending
 *			upon whether a particular type of information has
 *			been requested or not.  Only those types of information
 *			which have been requested (via the `I <info>' command
 *			at the user interface) will have their column headings
 *			printed.  The column heading is simply the name of the
 *			data as given by `info_tok[i]'.
 *
 *****************************************************************************/
void infotitles()
{
  short i, j;
  boolean firsti;
  firsti=TRUE;
  for (i=0; i<MINFO; i++) {
    if (info_list[i]) {
      if (!firsti) fprintf(ff,"\t");
      if (strcmp(info_tok[i],"rhoside")==0) {
        for (j=0; j<20; j++) fprintf(ff,"rho(%d)\t",j);
        fprintf(ff,"rho(20)");
      }
      else if (strcmp(info_tok[i],"bside")==0) {
        for (j=3; j<=20; j++) fprintf(ff,"%d-border\t",j);
        fprintf(ff,">20-border");
      }
      else if (strcmp(info_tok[i],"arean")==0) {
        for (j=0; j<20; j++) fprintf(ff,"area%d\t",j);
        fprintf(ff,"area20");
      }
      else {
        fprintf(ff,"%s",info_tok[i]);
      }
      firsti=FALSE;
    }
  }
  fprintf(ff,"\n");
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     I N F O W R I T E ( )  *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	infowrite()
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		This subroutine is the main engine for printing out
 *			information or statistics about the network.  The
 *			information is printed out in tab delimited format to
 *			the file "info.pb" (which must previously be opened by
 *			the routine `infoopen() ).  Only that information is
 *			printed out for which `info_list[i]' is equal to TRUE.
 *			This gives the user a mechanism for selecting which
 *			statistical information is to be recorded.  The type
 *			of information is identified by `info_tok[i]'.  A
 *			simple series of comparisons with `info_tok[i]'
 *			determines what to print out.  Some information
 *			requires a simple calculation to be performed, and
 *			there are a handful of subroutines below which
 *			calculate this data.
 *
 *****************************************************************************/
void infowrite()
{
  short i, j;
  short rho[21];
  REAL mu, arean[21], f, nbar;
  boolean firsti, wflag, fflag;
  REAL muarea(), phifn(), nbarfn(), efracfn(), zfn(), arootfn();
  void calcrho(), calcrhob(), calcarean();
  firsti=TRUE;
  for (i=0; i<MINFO; i++) {
    if (info_list[i]) {
      if (!firsti) fprintf(ff,"\t");
      wflag=info_wlist[i] & TOSCREEN;
      fflag=info_wlist[i] & TOFILE;
      if (strcmp(info_tok[i],"muside")==0) {
        calcrho(rho);
        nbar=nbarfn(rho);
        /* Remember to include bubbles in this calculation! */
        mu=0.0;
        for (j=0; j<21; j++) mu += (REAL) (j-nbar)*(j-nbar)*rho[j];
        mu /= (REAL) nc+nbub;
        if (fflag) fprintf(ff,"%f",mu);
        if (wflag) printf("muside = %f\n",mu);
      }
      else if (strcmp(info_tok[i],"rhoside")==0) {
        calcrho(rho);
        if (fflag) {
          for (j=0; j<20; j++) fprintf(ff,"%d\t",rho[j]);
          fprintf(ff,"%d",rho[20]);
        }
        if (wflag) {
          printf("distribution of sides:\n");
          for (j=0; j<21; j++) printf("rho(%d) = %d\n",j,rho[j]);
        }
      }
      else if (strcmp(info_tok[i],"muarea")==0) {
        mu=muarea();
        if (fflag) fprintf(ff,"%f",mu);
        if (wflag) printf("muarea = %f\n",mu);
      }
      else if (strcmp(info_tok[i],"arean")==0) {
        calcarean(arean);
        if (fflag) {
          for (j=0; j<20; j++) fprintf(ff,"%f\t",arean[j]);
          fprintf(ff,"%f",arean[20]);
        }
        if (wflag) {
          printf("normalized area of n-cells:\n");
          for (j=0; j<21; j++) printf("area%d = %f\n",j,arean[j]);
        }
      }
      else if (strcmp(info_tok[i],"bside")==0) {
        calcrhob(rho);
        if (fflag) {
          for (j=0; j<18; j++) fprintf(ff,"%d\t",rho[j]);
          fprintf(ff,"%d",rho[18]);
        }
        if (wflag) {
          printf("distribution of border sides:\n");
          for (j=0; j<18; j++) printf("rho(%d) = %d\n",j+3,rho[j]);
          printf("rho(>20) = %d\n",rho[18]);
        }
      }
      else if (strcmp(info_tok[i],"netenergy")==0) {
        if (fflag) fprintf(ff,"%f",netenergy);
        if (wflag) printf("netenergy = %f\n",netenergy);
      }
      else if (strcmp(info_tok[i],"henckyeps")==0) {
        if (fflag) fprintf(ff,"%f",henckyeps);
        if (wflag) printf("henckyeps = %f\n",henckyeps);
      }
      else if (strcmp(info_tok[i],"elosscount")==0) {
        if (fflag) fprintf(ff,"%d",elosscount);
        if (wflag) printf("elosscount = %d\n",elosscount);
      }
      else if (strcmp(info_tok[i],"bpinchcount")==0) {
        if (fflag) fprintf(ff,"%d",bpinchcount);
        if (wflag) printf("bpinchcount = %d\n",bpinchcount);
      }
      else if (strcmp(info_tok[i],"phi")==0) {
        f=phifn();
        if (fflag) fprintf(ff,"%f",f);
        if (wflag) printf("phi = %f\n",f);
      }
      else if (strcmp(info_tok[i],"nbar")==0) {
        calcrho(rho);
        nbar=nbarfn(rho);
        if (fflag) fprintf(ff,"%f",nbar);
        if (wflag) printf("nbar = %f\n",nbar);
      }
      else if (strcmp(info_tok[i],"edgefrac")==0) {
        f=efracfn();
        if (fflag) fprintf(ff,"%f",f);
        if (wflag) printf("edgefrac = %f\n",f);
      }
      else if (strcmp(info_tok[i],"z")==0) {
        calcrho(rho);
        f=zfn(rho);
        if (fflag) fprintf(ff,"%f",f);
        if (wflag) printf("z = %f\n",f);
      }
      else if (strcmp(info_tok[i],"arootbar")==0) {
        f=arootfn();
        if (fflag) fprintf(ff,"%f",f);
        if (wflag) printf("arootbar = %f\n",f);
      }
      firsti=FALSE;
    }
  }
  fprintf(ff,"\n");
  fflush(ff);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     * C A L C R H O ( )    *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	calcrho(short rho[] )
 *
 *	Arguments:	rho[n]	= the number of cells which have `n' sides
 *
 *	Return value:	none
 *
 *	Action:		Calculates the distribution of sides of cells (not
 *			normalised).  After a call to `calcrho(rho)', the
 *			array `rho[n]' contains the number of cells which
 *			have `n' sides.  At present the maximum number of
 *			sides considered is 20.
 *
 *****************************************************************************/
void calcrho(rho)
short rho[];
{
  short j, jj;
  for (j=0; j<21; j++) rho[j]=0;
  for (jj=0; jj<nc; jj++) {
    j=clist[jj];
    rho[ncsides[j]]++;
  }
  rho[0]=nbub;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  N B A R F N ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL nbarfn(short rho[])
 *
 *	Arguments:	rho[]	= value `rho[n]' is the number of cells with
 *				`n' sides
 *
 *	Return value:	average number of sides of a cell
 *
 *	Action:		Returns the average number of sides of a cell, given
 *			the distribution of cell sides `rho[n]' as an
 *			argument.  In a dry froth network this number will
 *			always be exactly `nbar = 6'.  However as the froth
 *			gets wetter (and 4-sided borders become possible) the
 *			value of `nbar' falls to a lower value.
 *
 *****************************************************************************/
REAL nbarfn(rho)
short rho[];
{
  short j;
  REAL nbar;
  nbar=0.0;
  for (j=0; j<21; j++) nbar += (REAL) j*rho[j];
  nbar /= (REAL) nc;
  return nbar;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *        Z F N ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL zfn(short rho[])
 *
 *	Arguments:	rho[]	= value `rho[n]' is the number of cells with
 *				`n' sides
 *
 *	Return value:	coordination number `z' of cells integral to network
 *
 *	Action:		Returns the coordination number of cells in the
 *			network.  This is almost identical to what is
 *			calculated by `nbarfn()' (the average number of cell
 *			sides) except that in this case we do not count zero
 *			or one-sided cells.  Therefore this statistic is more
 *			likely to converge to the value 4 as the froth
 *			approaches break up.
 *
 *****************************************************************************/
REAL zfn(rho)
short rho[];
{
  /* This fn returns the average number of sides to a cell (co-ordination) */
  /* *not* including zero or one sided cells */
  short j;
  REAL z;
  z=0.0;
  for (j=2; j<21; j++) z += (REAL) j*rho[j];
  z /= (REAL) (nc-(rho[0]+rho[1]));
  return(z);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  M U A R E A ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL muarea()
 *
 *	Arguments:	none
 *
 *	Return value:	Returns the second moment of the distribution of cell
 *			areas.
 *
 *****************************************************************************/
REAL muarea()
{
  short i, ii;
  REAL abar, sum, sqf();
  abar=boxwid*boxhgt/((REAL) nc);
  sum=0.0;
  for (ii=0; ii<nc; ii++) {
    i=clist[ii];
    sum += sqf(carea[i]+darea[i]-abar);
  }
  return( sum/(abar*abar*nc) );
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     C A L C A R E A N ( )  *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	calcarean(REAL arean[])
 *
 *	Arguments:	arean[]	= distribution of cell areas.  The value
 *				given by `arean[n]' is the average area of
 *				the n-sided cells divided by the area per cell
 *				`abar'.
 *				(Note that `abar' = network area/no. of cells).
 *
 *	Return value:	none
 *
 *	Action:		Calculates the area distribution `arean[n]'.
 *
 *****************************************************************************/
void calcarean(arean)
REAL arean[21];
{
  short i, ii, rho[21];
  REAL abar;
  void calcrho();
  abar=boxwid*boxhgt/((REAL) nc);
  /* First get the distribution of cell sides, `rho[n]' */
  calcrho(rho);
  for (i=0; i<=20; i++) arean[i]=0;
  for (ii=0; ii<nc; ii++) {
    i=clist[ii];
    arean[ncsides[i]] += carea[i]+darea[i];
  }
  for (i=0; i<=20; i++)
    if (rho[i]>0) arean[i] /= abar*rho[i];
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *C A L C R H O B ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	calcrhob(short rho[])
 *
 *	Arguments:	rho[]	= distribution of the sides of Plateau borders
 *
 *	Return value:	none
 *
 *	Action:		Returns the distribution of the number of sides of
 *			Plateau borders, where `rho[n-3]' records the number
 *			of Plateau borders wigh `n' edges.  The total of all
 *			monster sized Plateau bordes (with number of sides
 *			`n' > 20) goes into `rho[18]'.
 *
 *****************************************************************************/
void calcrhob(rho)
short rho[];
{
  short j, jj, ns;
  for (j=0; j<=20; j++) rho[j]=0;
  for (jj=0; jj<nb; jj++) {
    j=blist[jj];
    if ((ns=nbsides[j]) <= 20) rho[ns-3]++;
    else rho[18]++;
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     * E F R A C F N ( )    *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL efracfn()
 *
 *	Arguments:	none
 *
 *	Return value:	Using the concept of the decoration lemma, this
 *			subroutine estimates the average fraction of a
 *			cell-cell edge which lies buried inside its adjacent
 *			Plateau border (value lies in the range 0 to 1).
 *			Can be a useful alternative measure of wetness for
 *			foams which are not too wet.
 *
 *****************************************************************************/
REAL efracfn()
{
  short i, ii, validnv, b;
  REAL efrac, x1, y1, x2, y2, x11, y11, x12, y12, r, p1, p2, xp, yp, l, ll;
  REAL carclen();
  void vnbrxy();
  efrac=0.0; validnv=0;
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    if (nbsides[b=cadj[i][0]]==3) {
      validnv++;
      x1=vx[i]; y1=vy[i]; vnbrxy(i,0,&x2,&y2);
      vnbrxy(i,1,&x11,&y11); vnbrxy(i,2,&x12,&y12);
      p1=cp[cadj[i][1]]; p2=cp[cadj[i][2]];
      /* This is a crude way of finding the centre of a PB */
      xp=(x1+x11+x12)/3.0; yp=(y1+y11+y12)/3.0;
      l=0.5*carclen(x1,y1,x2,y2,p1,p2);
      ll=carclen(xp,yp,x1,y1,p1,p2);
      efrac += ll/(l+ll);
    }
  }
  efrac /= (REAL) validnv;
  return(efrac);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     * A R O O T F N ( )    *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	REAL arootfn()
 *
 *	Arguments:	none
 *
 *	Return value:	Returns the average of the square root of each cell
 *			area, i.e. < sqrt(A_cell) > which is a useful
 *			characteristic linear dimension of the network.
 *
 *****************************************************************************/
REAL arootfn()
{
  short i, ii, c;
  REAL ra, cellarea();
  boolean found[MCELL], arcbrk;
  for (i=0; i<onc; i++) found[i]=FALSE;
  ra=0.0;
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    if (!found[c=cadj[i][1]]) {
      found[c]=TRUE;
      ra += sqrt(cellarea(i,1,&arcbrk));
    }
  }
  ra /= (REAL) nc;
  return(ra);
}
 
/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     I N F O C L O S E ( )  *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	infoclose()
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		Speaks for itself...complementary to `infoopen()'.
 *
 *****************************************************************************/
void infoclose()
{
  fclose(ff);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  S E T I N F O ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	setinfo(char *tok, boolean wflag, fflag)
 *
 *	Arguments:	tok[]	= name of piece of information requested
 *			wflag	= TRUE if you want to write to the screen
 *			fflag	= TRUE if you want to write to file "info.pb"
 *
 *	Return value:	none
 *
 *	Action:		Requests that the information type named by `tok[]'
 *			be included in the data which is written by the routine
 *			`infowrite()'.  The arguments `wflag' and `fflag'
 *			further qualify this by saying if you want the info
 *			written to `stdout' and/or to the file `info.pb',
 *			respectively.  This is used to implement the command
 *			`I <infotype>' which is part of the command interface.
 *			Note that this routine should not be called while the
 *			file `info.pb' is open or else the number of columns
 *			in the file would change half way through, making it
 *			difficult to read.
 *
 *****************************************************************************/
void setinfo(tok,wflag,fflag)
char *tok;
boolean wflag, fflag;
{
  void plerror();
  short i, j;
  for (i=0; i<MINFO; i++) {
    if (strcmp(tok,info_tok[i])==0) {
      info_list[i]=TRUE;
      info_wlist[i] |= (wflag) ? TOSCREEN : 0;
      info_wlist[i] |= (fflag) ? TOFILE : 0;
      return;
    }
  }
  plerror("requested unrecognized information");
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *      C A D J M ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	cadjm()
 *
 *	Arguments:	 
 *
 *	Return value:	none
 *
 *	Action:		Traveses all cells and for each cell, walks it's	
 *			perimeter, recording all neighbours in an 
 *			adjacency matrix.
 *			
 *			This matrix is then checked for symmetry and
 *			then output into a sequentially numbered file.
 *			
 *
 *****************************************************************************/
void cadjm()
{
  FILE * matout; // file pointer to write matrix out
  static int callnumb = 0; // counter to give each file a unique name
  char filename[20]; // string for filename
  sprintf(filename,"mat%03d.dat",callnumb);
  if( (matout = fopen(filename,"w")) == NULL){
    plerror("error opening file for adj mat output");
    return;
  }
  /*
   * Allocate the adjacency matrix and set up simple 2 dimensional array
   */
  int ** adjmat = calloc(nc,sizeof(int*));
  int * data = calloc(nc*nc, sizeof(int));
  int c,i,ii,k;
  for(c=0;c<nc;c++) adjmat[c] = data + c*nc;

  for(c=0;c<nc;c++){ // for each cell
    adjmat[c][c] = 1;
    for (ii=0; ii<nv; ii++) { // find a vertex on that cell
      i=vlist[ii];
      if (c==cadj[i][k=1]) break;
      if (c==cadj[i][k=2]) break;
    }
    if (k==2) { i=vnbr[i][0]; k=1; } // chose orientation
    ii=i; // remember starting point
    do {
      k = cadj[i][2]; // cell on 'left' is a neighbour
      adjmat[c][k] = 1;
      i=vnbr[i][2]; i=vnbr[i][0]; // walk perimeter
    } while (i!=ii); // until you return to start
  }

  /*
   * check that resulting matrix is symetric
   */
  for(c=0;c<nc;c++){
    for(k=c;k<nc;k++){
      if(adjmat[c][k] != adjmat[k][c])
        plerror("adjacency matrix is not symmetric!");
    }
  }

  for(c=0;c<nc;c++){
    for(k=0;k<nc;k++){
      fprintf(matout,"%d\t",adjmat[c][k]);
    }
    fprintf(matout,"\n");
  }
  /* 
   *  free all allocated memory and close all opened filepointers
   */
  free(data);
  free(adjmat);
  fclose(matout);
  callnumb++;
}
