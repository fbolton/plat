#include "include.h"
#ifdef IRIS
#include <device.h>
#endif

/*********************************/
/* 'mgl' - 'my graphics library' */
/*********************************/

#ifdef IRIS
static REAL irxo= IRISXO, iryo= IRISYO, irxs= IRISXS, irys= IRISYS;
long preflx= IRISPREFLX, prefrx= IRISPREFRX, prefby= IRISPREFBY,
     prefty= IRISPREFTY;
#endif
#if defined(HPGL) || defined(ORIHPGL)
static REAL hpglxo= HPGLXO, hpglyo= HPGLYO, hpglxs= HPGLXS, hpglys= HPGLYS;
#endif
#ifdef TEK
static REAL tekxo= TEKXO, tekyo= TEKYO, tekxs= TEKXS, tekys= TEKYS;
#endif
#ifdef POSTSCRIPT
static REAL postxo= POSTXO, postyo= POSTYO, postxs= POSTXS, postys= POSTYS;
#endif
#ifdef X_WINDOWS
static REAL x_win_xo= X_WIN_XO, x_win_yo= X_WIN_YO,
            x_win_xs= X_WIN_XS, x_win_ys= X_WIN_YS;
static Widget wtop;
#endif

#define INTOF(a) ((int) (0.5+(a)) )

/* This line is included for the HP-GL Graphics Language. */
/* (with the option of specifying the O'Reilly plotter)   */
#if defined(HPGL) || defined(ORIHPGL)
FILE *hpgl;
#endif

/* These lines are included to generate a Tektronix file. */
#ifdef TEK
FILE *tek;
short tekcount=0;
#endif

/* This line is included to generate a postscript file. */
#ifdef POSTSCRIPT
FILE *post;
#endif

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *      Z O O M ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	zoom()
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		When IRIS graphics are enabled this subroutine is
 *			used to implement a `zoom' feature.  This subroutine
 *			reads queued mouse events (using IRIS's own built-in
 *			subroutines) and changes the scale of the IRIS
 *			graphic window accordingly.  A LEFTMOUSE button event
 *			causes the screen to zoom out, centred on the position
 *			of the mouse click.  A MIDDLEMOUSE button event causes
 *			the origin to be repositioned at the point where the
 *			mouse was pointing.  A RIGHTMOUSE button event causes
 *			the screen to zoom in, centred on the mouse position.
 *			The global variables `irxo, iryo, irxs, irys' are
 *			changed to reflect the rescaling.  Zooming does not
 *			take effect until the IRIS window is replotted with
 *			these new parameters.
 *
 *			Note that the event queue for these mouse events is
 *			set up by the routine `mglreset()' whenever IRIS
 *			graphics are enabled.
 *
 *****************************************************************************/
#ifdef IRIS
void zoom()
{
  short mr, mx, my;
  REAL dxo, dyo;
  Device dev;
/*Screencoord lx, rx, by, ty;*/
/*qreset();   flushes the input buffer */
  while (ginteract && (dev=qtest())) {
    if (dev==RIGHTMOUSE || dev==MIDDLEMOUSE || dev==LEFTMOUSE) {
      qread(&mr);
      qread(&mx); qread(&my);
      if (mr==1) {
        dxo= (((REAL)mx-(REAL)preflx)-irxo)/irxs;
        dyo= (((REAL)my-(REAL)prefby)-iryo)/irys;
        if (dev==RIGHTMOUSE) {
          irxs *= 1.5; irys *= 1.5;
        }
        else if (dev==LEFTMOUSE) {
          irxs /= 1.5; irys /= 1.5;
        }
        irxo= IRISXO-irxs*dxo; iryo= IRISYO-irys*dyo;
      }
    }
    else {
      qread(&mr);
    }
  }
}
#endif

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *      Z O O M ( )     *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	zoom()
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		Another implementaton of the `zoom()' function, this
 *			time for the X windows interface (many of the details
 *			of the X windows interface are contained in the file
 *			`xplat.c').
 *
 *			The effect of the routine is to rescale the global
 *			parameters `x_win_xo, x_win_yo, x_win_xs, x_win_ys' so
 *			that the next time the foam is redrawn within the X
 *			window it will have been rescaled.  As with the other
 *			version of the `zoom()' function, a left button zooms
 *			out, a middle button repositions the origin and a right
 *			button press zooms in.
 *
 *****************************************************************************/
#ifdef X_WINDOWS
void zoom(int button, int x, int y)
{
  static REAL zoom_factor = 1.1;
  REAL dxo, dyo;

  dxo = (((REAL) x) - x_win_xo)/x_win_xs;
  dyo = (((REAL) y) - x_win_yo)/x_win_ys;
  switch (button) {
    case 1 : x_win_xs *= zoom_factor;
	x_win_ys *= zoom_factor;
	break;
    case 2 : break;
    case 3 : x_win_xs /= zoom_factor;
	x_win_ys /= zoom_factor;
	break;
  }
  x_win_xo = X_WIN_XO - x_win_xs*dxo;
  x_win_yo = X_WIN_YO - x_win_ys*dyo;
#ifdef HPGL
  switch (button) {
    case 1 : hpglxs *= zoom_factor;
	hpglys *= zoom_factor;
	break;
    case 2 : break;
    case 3 : hpglxs /= zoom_factor;
	hpglys /= zoom_factor;
	break;
  }
  hpglxo = HPGLXO - hpglxs*dxo;
  hpglyo = HPGLYO - hpglys*dyo;
#endif
#ifdef TEK
  switch (button) {
    case 1 : tekxs *= zoom_factor;
	tekys *= zoom_factor;
	break;
    case 2 : break;
    case 3 : tekxs /= zoom_factor;
	tekys /= zoom_factor;
	break;
  }
  tekxo = TEKXO - tekxs*dxo;
  tekyo = TEKYO - tekys*dyo;
#endif
#ifdef POSTSCRIPT
  switch (button) {
    case 1 : postxs *= zoom_factor;
	postys *= zoom_factor;
	break;
    case 2 : break;
    case 3 : postxs /= zoom_factor;
	postys /= zoom_factor;
	break;
  }
  postxo = POSTXO - postxs*dxo;
  postyo = POSTYO - postys*dyo;
#endif
}
#endif

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  M G L R E S E T ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	mglreset()
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		Resets all scaling variables and the position of the
 *			origin so that the picture of the foam will be
 *			displayed/plotted/etc. in its usual position at the
 *			usual size.
 *
 *			In the case where IRIS graphics are used, it also
 *			initialises the IRIS event queue which is needed to
 *			detect mouse button events on the IRIS machine.
 *
 *****************************************************************************/
void mglreset()
{
#ifdef IRIS
  irxo= IRISXO; iryo= IRISYO; irxs= IRISXS; irys= IRISYS;
  if (ginteract) {
    preflx= IRISPREFLX; prefrx= IRISPREFRX;
    prefby= IRISPREFBY; prefty= IRISPREFTY;
    prefposition(preflx,prefrx,prefby,prefty);
#if defined(DEBUG) && defined(IRIS)
    foreground();
#endif
    winopen("foam");
    doublebuffer();
    gconfig();
    color(BLACK);
    clear();
    color(CYAN);
    swapbuffers();
    qdevice(RIGHTMOUSE);
    tie(RIGHTMOUSE,MOUSEX,MOUSEY);
    qdevice(LEFTMOUSE);
    tie(LEFTMOUSE,MOUSEX,MOUSEY);
    qdevice(MIDDLEMOUSE);
    tie(MIDDLEMOUSE,MOUSEX,MOUSEY);
  }
#endif
#if defined(HPGL) || defined(ORIHPGL)
  hpglxo= HPGLXO; hpglyo= HPGLYO; hpglxs= HPGLXS; hpglys= HPGLYS;
#endif
#ifdef TEK
  tekxo= TEKXO; tekyo= TEKYO; tekxs= TEKXS; tekys= TEKYS;
#endif
#ifdef POSTSCRIPT
  postxo= POSTXO; postyo= POSTYO; postxs= POSTXS; postys= POSTYS;
#endif
#ifdef X_WINDOWS
  wtop = bitmap.frothPicture;
  x_win_xo= X_WIN_XO; x_win_yo= X_WIN_YO;
  x_win_xs= X_WIN_XS; x_win_ys= X_WIN_YS;
#endif
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *    M G L I N I T ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	mglinit()
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		The action of this subroutine varies depending upon
 *			the graphics package selected.  Broadly speaking it
 * 			issues whatever initialisation commands are needed
 *			by the graphics package before a new picture is
 *			drawn.  Generally, the structure of any drawing
 *			sequence would be:
 *				mglinit();
 *				(miscellaneous graphics commands...);
 *				mglclose();
 *			So that `mglinit()' and `mglclose()' are used as a 
 *			matching pair of commands.
 *
 *****************************************************************************/
void mglinit()
{
#ifdef IRIS
  if (ginteract) {
    color(BLACK);
    clear();
    color(CYAN);
  }
#endif
#ifdef HPGL
  if (mgl_hpgl_flag) {
    hpgl=fopen(mgl_hpgl_filename,"w");
    fprintf(hpgl,"IN;SP1;\n");
    fprintf(hpgl,"IP%.0f,%.0f,%.0f,%.0f;SC%.0f,%.0f,%.0f,%.0f;\n",hpglxo,hpglyo,hpglxo+hpglxs,hpglyo+hpglys,0.0,1.0,0.0,1.0);
  }
#endif
#ifdef TEK
  if (mgl_tek_flag) {
    tekcount=0;
    tek=fopen(mgl_tek_filename,"w");
    fprintf(tek,"%c%c",0x1b,0x0c);
    tekcount += 2;
  }
#endif
#ifdef POSTSCRIPT
  if (mgl_ps_flag) {
    post=fopen(mgl_ps_filename,"w");
    fprintf(post,"%%!PS-Adobe-2.0\n");
    fprintf(post,"%%%%Title: foam\n");
    fprintf(post,"%%%%EndComments\n");
    fprintf(post,"/m /moveto load def\n");
    fprintf(post,"/l /lineto load def\n");
    fprintf(post,"%%%% define standard linewidth\n");
    fprintf(post,"\n");
    fprintf(post,"0.25 0.25 scale\n");
    fprintf(post,"%.1f setlinewidth\n\n",POSTLINEWID);
  }
#endif
#ifdef X_WINDOWS
  wtop = bitmap.frothPicture;
  XFillRectangle(
	XtDisplay(wtop),
	bitmap.froth_bitmap,
	bitmap.undraw_gc,
	0, 0,
	bitmap.pixmap_width+1,
	bitmap.pixmap_height+1);
#endif
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  M G L P O I N T ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	mglpoint(REAL x1, REAL y1)
 *
 *	Arguments:	(x1, y1) = Coordinate of point
 *
 *	Return value:	none
 *
 *	Action:		Plots a point (on all of the current graphics devices.)(remains to be seen)
 *			For plotting parameters see file `include.h'.
 *
 *****************************************************************************/
void mglpoint(x1,y1)
REAL x1,y1;
{
/*
 * Commented out everything copied from mglline 
 * that is not postscript or xwindows
 *
 * #ifdef IRIS
  REAL xy[2];
#endif
#ifdef TEK
  short ix1, iy1, ix2, iy2;
  char b1, b2, b3, b4;
#endif
#if defined(HPGL) || defined(ORIHPGL)
  short i;
  REAL dx, dy;
#endif
#ifdef IRIS
  if (ginteract) {
    bgnline();
    xy[0]=irxo+irxs*x1; xy[1]=iryo+irys*y1;
#ifdef DOUBLEPRECISION
    v2d(xy);
#else
    v2f(xy);
#endif
    xy[0]=irxo+irxs*x2; xy[1]=iryo+irys*y2;
#ifdef DOUBLEPRECISION
    v2d(xy);
#else
    v2f(xy);
#endif
    endline();
  }
#endif
#if defined(HPGL)
  if (mgl_hpgl_flag) {
[>
    x1=hpglxo+hpglxs*x1; y1=hpglyo+hpglys*y1;
    x2=hpglxo+hpglxs*x2; y2=hpglyo+hpglys*y2;
    fprintf(hpgl,"PU%.0f,%.0f;",x1,y1);
    dx=(x2-x1)/20.0; dy=(y2-y1)/20.0;
    for (i=1; i<=20; i++) {
      fprintf(hpgl,"PD%.0f,%.0f;",x1+i*dx,y1+i*dy);
    }
<]
    fprintf(hpgl,"PU%.0f,%.0f;PD%.0f,%.0f;\n",x1,y1,x2,y2);
  }
#endif
#ifdef TEK
  [>
   *  Lots of Black Art here -- I admit it's incomprehensible,
   *  and I don't know what it does myself any more!
   <]
  if (mgl_tek_flag) {
    ix1= (short) (tekxo+tekxs*x1); iy1= (short) (tekyo+tekys*y1);
    ix2= (short) (tekxo+tekxs*x2); iy2= (short) (tekyo+tekys*y2);
    if (ix1>=0 && ix1<1024 && iy1>=0 && iy1<1024
     && ix2>=0 && ix2<1024 && iy2>=0 && iy2<1024) {
      b1=((char) ((iy1&0x3e0)>>5)|0x20);
      b2=((char) (iy1&0x1f)|0x60);
      b3=((char) ((ix1&0x3e0)>>5)|0x20);
      b4=((char) (ix1&0x1f)|0x40);
      fprintf(tek,"%c%c%c%c%c",0x1d,b1,b2,b3,b4);
      b1=((char) ((iy2&0x3e0)>>5)|0x20);
      b2=((char) (iy2&0x1f)|0x60);
      b3=((char) ((ix2&0x3e0)>>5)|0x20);
      b4=((char) (ix2&0x1f)|0x40);
      [>
       *  The following is a bizarre kludge which proved to
       *  be necessary when this code was first developed.
       *  A whole brace of `no-ops' are sent to the Tektronix
       *  display to give the terminal some time to think!
       <]
      fprintf(tek,"%c%c%c%c%c%c%c%c%c%c%c%c",
        b1,b2,b3,b4,0x16,0x16,0x16,0x16,0x16,0x16,0x16,0x16);
      if ((tekcount += 17)>490) { fprintf(tek,"%c",0x0a); tekcount=0; }
    }
  }
#endif*/
#ifdef POSTSCRIPT
  if (mgl_ps_flag) {
    fprintf(post,"%.1f %.1f 5 0 360 newpath arc\nfill\n",postxo+postxs*x1,postyo+postys*y1);
    /*fprintf(post,"%.1f %.1f l\n",postxo+postxs*x2,postyo+postys*y2);*/
  }
#endif
#ifdef X_WINDOWS
  XDrawPoint(XtDisplay(wtop),bitmap.froth_bitmap,bitmap.draw_gc,
	    INTOF(x_win_xo+x_win_xs*x1), INTOF(x_win_yo+x_win_ys*y1) );
#endif  
}


/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *    M G L L I N E ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	mglline(REAL x1, REAL y1, REAL x2, REAL y2)
 *
 *	Arguments:	(x1, y1) = Coordinate of first endpoint of line
 *			(x2, y2) = Coordinate of second endpoin of line
 *
 *	Return value:	none
 *
 *	Action:		Plots a line on all of the current graphics devices.
 *			For plotting parameters see file `include.h'.
 *
 *****************************************************************************/
void mglline(x1,y1,x2,y2)
REAL x1,y1,x2,y2;
{
#ifdef IRIS
  REAL xy[2];
#endif
#ifdef TEK
  short ix1, iy1, ix2, iy2;
  char b1, b2, b3, b4;
#endif
#if defined(HPGL) || defined(ORIHPGL)
  short i;
  REAL dx, dy;
#endif
#ifdef IRIS
  if (ginteract) {
    bgnline();
    xy[0]=irxo+irxs*x1; xy[1]=iryo+irys*y1;
#ifdef DOUBLEPRECISION
    v2d(xy);
#else
    v2f(xy);
#endif
    xy[0]=irxo+irxs*x2; xy[1]=iryo+irys*y2;
#ifdef DOUBLEPRECISION
    v2d(xy);
#else
    v2f(xy);
#endif
    endline();
  }
#endif
#if defined(HPGL)
  if (mgl_hpgl_flag) {
/*
    x1=hpglxo+hpglxs*x1; y1=hpglyo+hpglys*y1;
    x2=hpglxo+hpglxs*x2; y2=hpglyo+hpglys*y2;
    fprintf(hpgl,"PU%.0f,%.0f;",x1,y1);
    dx=(x2-x1)/20.0; dy=(y2-y1)/20.0;
    for (i=1; i<=20; i++) {
      fprintf(hpgl,"PD%.0f,%.0f;",x1+i*dx,y1+i*dy);
    }
*/
    fprintf(hpgl,"PU%.0f,%.0f;PD%.0f,%.0f;\n",x1,y1,x2,y2);
  }
#endif
#ifdef TEK
  /*
   *  Lots of Black Art here -- I admit it's incomprehensible,
   *  and I don't know what it does myself any more!
   */
  if (mgl_tek_flag) {
    ix1= (short) (tekxo+tekxs*x1); iy1= (short) (tekyo+tekys*y1);
    ix2= (short) (tekxo+tekxs*x2); iy2= (short) (tekyo+tekys*y2);
    if (ix1>=0 && ix1<1024 && iy1>=0 && iy1<1024
     && ix2>=0 && ix2<1024 && iy2>=0 && iy2<1024) {
      b1=((char) ((iy1&0x3e0)>>5)|0x20);
      b2=((char) (iy1&0x1f)|0x60);
      b3=((char) ((ix1&0x3e0)>>5)|0x20);
      b4=((char) (ix1&0x1f)|0x40);
      fprintf(tek,"%c%c%c%c%c",0x1d,b1,b2,b3,b4);
      b1=((char) ((iy2&0x3e0)>>5)|0x20);
      b2=((char) (iy2&0x1f)|0x60);
      b3=((char) ((ix2&0x3e0)>>5)|0x20);
      b4=((char) (ix2&0x1f)|0x40);
      /*
       *  The following is a bizarre kludge which proved to
       *  be necessary when this code was first developed.
       *  A whole brace of `no-ops' are sent to the Tektronix
       *  display to give the terminal some time to think!
       */
      fprintf(tek,"%c%c%c%c%c%c%c%c%c%c%c%c",
        b1,b2,b3,b4,0x16,0x16,0x16,0x16,0x16,0x16,0x16,0x16);
      if ((tekcount += 17)>490) { fprintf(tek,"%c",0x0a); tekcount=0; }
    }
  }
#endif
#ifdef POSTSCRIPT
  if (mgl_ps_flag) {
    fprintf(post,"%.1f %.1f m\n",postxo+postxs*x1,postyo+postys*y1);
    fprintf(post,"%.1f %.1f l\n",postxo+postxs*x2,postyo+postys*y2);
  }
#endif
#ifdef X_WINDOWS
  XDrawLine(XtDisplay(wtop),bitmap.froth_bitmap,bitmap.draw_gc,
	    INTOF(x_win_xo+x_win_xs*x1), INTOF(x_win_yo+x_win_ys*y1),
	    INTOF(x_win_xo+x_win_xs*x2), INTOF(x_win_yo+x_win_ys*y2) );
#endif  
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *    M G L A R C ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	mglarc(REAL xc, REAL yc, REAL r, REAL th1, REAL th2)
 *
 *	Arguments:	(xc, yc) = Centre coordinates of the arc
 *			r	= Radius of arc
 *			th1	= Initial angle of arc (anti-clockwise)
 *			th2	= Final angle of arc (anti-clockwise)
 *
 *	Return value:	none
 *
 *	Action:		Draws an arc centred at (xc, yc) of radius `r'.  The
 *			arc is drawn between the angles `th1' and `th2' moving
 *			in an anti-clockwise direction.
 *			For plotting parameters see file `include.h'.
 *			
 *
 *****************************************************************************/
void mglarc(xc, yc, r, th1, th2)
REAL xc, yc, r, th1, th2;
{
#ifdef TEK
  REAL dtheta;
  short i, n, ix, iy, ix1, iy1, ix2, iy2;
  char b1, b2, b3, b4;
  REAL x1,y1,x2,y2;
#endif
#if defined(POSTSCRIPT)
  REAL dtheta;
  short i, n;
  REAL x1,y1,x2,y2;
#endif
#ifdef IRIS
  if (ginteract) {
#ifdef DOUBLEPRECISION
    arc((float) irxo+irxs*xc,(float) iryo+irys*yc,(float) irxs*r,
      (float) 1800.0*th1/PI,(float) 1800.0*th2/PI);
#else
    arc(irxo+irxs*xc, iryo+irys*yc, irxs*r, 1800.0*th1/PI, 1800.0*th2/PI);
#endif
  }
#endif
#if defined(HPGL)
  if (mgl_hpgl_flag) {
    fprintf(hpgl,"PU%.0f,%.0f;PD;",(xc+r*cos(th1)), (yc+r*sin(th1)));
    fprintf(hpgl,"AA%.0f,%.0f,%.1f;\n", xc, yc, 180.0*(th2-th1)/PI);
  }
#endif
#ifdef TEK
  if (mgl_tek_flag) {
    n= TEKARCNUM;
    dtheta=(th2-th1)/((REAL) n);
    x2= xc+r*cos(th1);
    y2= yc+r*sin(th1);
    for (i=1; i<=n; i++) {
      x1=x2; y1=y2;
      x2= xc+r*cos(th1+i*dtheta);
      y2= yc+r*sin(th1+i*dtheta);
      /* Now draw a line between x1,y1 and x2,y2 */
      ix1= (short) (tekxo+tekxs*x1); iy1= (short) (tekyo+tekys*y1);
      ix2= (short) (tekxo+tekxs*x2); iy2= (short) (tekyo+tekys*y2);
      if (ix1>=0 && ix1<1024 && iy1>=0 && iy1<1024
       && ix2>=0 && ix2<1024 && iy2>=0 && iy2<1024) {
        b1=((char) ((iy1&0x3e0)>>5)|0x20);
        b2=((char) (iy1&0x1f)|0x60);
        b3=((char) ((ix1&0x3e0)>>5)|0x20);
        b4=((char) (ix1&0x1f)|0x40);
        fprintf(tek,"%c%c%c%c%c",0x1d,b1,b2,b3,b4);
        b1=((char) ((iy2&0x3e0)>>5)|0x20);
        b2=((char) (iy2&0x1f)|0x60);
        b3=((char) ((ix2&0x3e0)>>5)|0x20);
        b4=((char) (ix2&0x1f)|0x40);
        fprintf(tek,"%c%c%c%c%c%c%c%c%c%c%c%c",
          b1,b2,b3,b4,0x16,0x16,0x16,0x16,0x16,0x16,0x16,0x16);
        if ((tekcount += 17)>490) { fprintf(tek,"%c",0x0a); tekcount=0; }
      }
    }
  }
#endif
#ifdef POSTSCRIPT
  if (mgl_ps_flag) {
    x1= xc+r*cos(th1);
    y1= yc+r*sin(th1);
    x2= xc+r*cos(th1);
    y2= yc+r*sin(th1);
    fprintf(post,"%.1f %.1f m\nstroke\n",postxo+postxs*x2,postyo+postys*y2);
    fprintf(post,"%.1f %.1f %.1f %.1f %.1f arc\n",postxo+postxs*xc,postyo+postys*yc,postxs*r,180.0*th1/PI,180.0*th2/PI);
  }
#endif
#ifdef X_WINDOWS
  XDrawArc(XtDisplay(wtop),bitmap.froth_bitmap,bitmap.draw_gc,
	    INTOF(x_win_xo+x_win_xs*(xc-r)), INTOF(x_win_yo+x_win_ys*(yc-r)),
	    INTOF(x_win_xs*2.0*r),           INTOF(x_win_ys*2.0*r),
	    (INTOF(-11520.0*th1/PI)),     (INTOF(-11520*(th2-th1)/PI)) );
#endif  
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *    M G L R E C T ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	mglrect(REAL x1, REAL y1, REAL x2, REAL y2)
 *
 *	Arguments:	(x1, y1) = Coordinates of one corner or rectangle
 *			(x2, y2) = Coordinates of diagonally opposite corner
 *
 *	Return value:	none
 *
 *	Action:		Plots a rectangle on all graphics devices where one
 *			corner is given by (x1, y1) and the diagonally
 *			opposite corner is given by (x2, y2).
 *			For plotting parameters see file `include.h'.
 *
 *****************************************************************************/
void mglrect(x1, y1, x2, y2)
REAL x1, y1, x2, y2;
{
#ifdef TEK
  void mglline();
#endif
#ifdef IRIS
  if (ginteract) {
#ifdef DOUBLEPRECISION
    rect((float) irxo+irxs*x1,(float) iryo+irys*y1,(float) irxo+irxs*x2,(float) iryo+irys*y2);
#else
    rect(irxo+irxs*x1, iryo+irys*y1, irxo+irxs*x2, iryo+irys*y2);
#endif
  }
#endif
#if defined(HPGL)
  if (mgl_hpgl_flag) {
    fprintf(hpgl,"PU%.0f,%.0f;", x1, y1);
    fprintf(hpgl,"EA%.0f,%.0f;", x2, y2);
  }
#endif
#ifdef TEK
  if (mgl_tek_flag) {
    mglline(x1,y1,x1,y2); mglline(x1,y2,x2,y2);
    mglline(x2,y2,x2,y1); mglline(x2,y1,x1,y1);
  }
#endif
#ifdef POSTSCRIPT
  if (mgl_ps_flag) {
    fprintf(post,"%.1f %.1f m\n%.1f %.1f l\n%.1f %.1f l\n%.1f %.1f l\n%.1f %.1f l\nclosepath\n",
          postxo+postxs*x1,postyo+postys*y1,postxo+postxs*x1,postyo+postys*y2,
          postxo+postxs*x2,postyo+postys*y2,postxo+postxs*x2,postyo+postys*y1,
          postxo+postxs*x1,postyo+postys*y1);
  }
#endif
#ifdef X_WINDOWS
  XDrawRectangle(XtDisplay(wtop),bitmap.froth_bitmap,bitmap.draw_gc,
	    (int) x_win_xo+x_win_xs*x1, (int) x_win_yo+x_win_ys*y1,
	    (int) x_win_xs*(x2-x1), (int) x_win_ys*(y2-y1));
#endif  
}
   
/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  M G L C L O S E ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	mglclose()
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		Finishes off the process of drawing a picture.  In the
 *			case of interactive graphics, it is at this point that
 *			the image is mapped to the display.  Used in
 *			conjunction with `mglinit()', for example:
 *				mglinit();
 *				(miscellaneous graphics commands...);
 *				mglclose();
 *			So that `mglinit()' and `mglclose()' are used as a 
 *			matching pair of commands.
 *
 *****************************************************************************/
void mglclose()
{
#ifdef IRIS
  if (ginteract) { swapbuffers(); }
#endif
#if defined(HPGL)
  if (mgl_hpgl_flag) {
    fprintf(hpgl,"PU;SP0;\n");
    fclose(hpgl);
  }
#endif
#ifdef TEK
  if (mgl_tek_flag) {
    fprintf(tek,"%c%c",0x1f,0x0a);
    fclose(tek);
  }
#endif
#ifdef POSTSCRIPT
  if (mgl_ps_flag) {
    fprintf(post,"\nstroke\nshowpage\n");
    fclose(post);
  }
#endif
#ifdef X_WINDOWS
  if (DefaultDepthOfScreen(XtScreen(wtop)) == 1)
    XCopyArea(XtDisplay(wtop),
	bitmap.froth_bitmap,
	XtWindow(wtop),
	DefaultGCOfScreen(XtScreen(wtop)), 0, 0,
	bitmap.pixmap_width,
	bitmap.pixmap_height,
	0, 0);
  else
    XCopyPlane(XtDisplay(wtop),
	bitmap.froth_bitmap,
	XtWindow(wtop),
	DefaultGCOfScreen(XtScreen(wtop)), 0, 0,
	bitmap.pixmap_width,
	bitmap.pixmap_height,
	0, 0, 1);
#endif
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *    M G L F O A M ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	mglfoam(REAL xorigin, REAL yorigin, REAL scale,
 *			        boolean leftfil, botfil,
 *			        REAL theta,
 *			        boolean hatch, centroid)
 *
 *	Arguments:	(xorigin, yorigin) = The centre of the foam plot
 *				relative to the display origin.  These are
 *				given in the same units as unscaled coords.
 *				(i.e. the size of the display window is
 *				assumed to be [-0.5,0.5]x[-0.5,0.5] ).
 *				(ignored if `scale == 0.0')
 *			scale	= Scales the size of the plot (usually 1.0)
 *				(defaults to 1.0 if scale == 0.0)
 *			leftfil	= Flag which says whether to plot any arcs
 *				which straddle the left edge of the enclosing
 *				rectangle (ignored if `scale == 0.0')
 *			botfile	= Flag which says whether to plot any arcs
 *				which straddle the bottom edge of the enclosing
 *				rectangle (ignored if `scale == 0.0')
 *			theta	= slope of the hatching lines which fill the
 *				Plateau borders (ignored if `hatch == FALSE')
 *			hatch	= Flag to indicate whether hatching of the
 *				Plateau borders is desired.
 *			centroid	= Flag to indicate whether drawing of centroids
 *				of bubbles is desired.
 *
 *	Return value:	none
 *
 *	Action:		This subroutine plots the foam network, and allows a
 *			few optional parameters to be specified.  A box of
 *			size 1x1 is drawn around the finished network.  The
 *			reason for the extra parameters `scale, xorigin,
 *			yorigin, leftfil, botfil' is to facilitate the
 *			plotting of many periodic boxes (e.g. when this
 *			subroutine is called via the `polyplot()' routine).
 *			Otherwise this routine is usually called indirectly
 *			via the routine `foamplot()'.
 *
 *			This routine is called in two possible modes:
 *			`rescaled' (when called by `polyplot()') and
 *			`non-rescaled' (when called by `foamplot()').
 *			By setting `scale = 0' the plot is deemed to be
 *			non-rescaled and the first five arguments are all
 *			ignored.
 *
 *			If `scale' is non-zero then this subroutine is in
 *			rescaled mode which allows the position of the network
 *			to be shifted and its size to be rescaled.  The coords
 *			(xorigin, yorigin) determine where the centre of the
 *			network is shifted to and all network coordinates are
 *			then multiplied by scale (relative to the new centre).
 *			Use of `leftfil' and `botfil' avoids avoids replotting
 *			arcs along the boundary when many periodic boxes are
 *			being plotted side by side.
 *
 *			The flag `hatch' indicates whether you want the
 *			Plateau borders to be shaded by hatching.  If so, then
 *			the hatching lines will be tilted at an angle given by
 *			the argument `theta'.
 *
 *****************************************************************************/
void mglfoam(xorigin,yorigin,scale,leftfil,botfil,theta,hatch,centroid)
REAL xorigin, yorigin, scale, theta;
boolean leftfil, botfil, hatch, centroid;
{
  short i, ii, k, j, j1, c1, c2, b;
  REAL x1, y1, x2, y2, p1, p2, alpha, xc, yc, r, th, th1, temp;
  REAL barcangle(), carcangle(), linangle(), linlen(), twopimod();
  boolean found[MBORD], straight, arcbrk, rescaled, la, arccentre(), larc();
  void mglline(), mglarc(), mglrect(), trans(), bubinbox(), mglhatch();
  void plerror();
  rescaled= (scale!=0);
  if (rescaled) {
    xorigin *= boxwid;
    yorigin *= boxhgt;
  }
  for (ii=0; ii<nv; ii++) {
    i=vlist[ii];
    /* Case of a cell-arc */
    x1=vx[i]; y1=vy[i];
    j=vnbr[i][0]; j1=(vper[i][0] & PERMASK);
    c1=cadj[i][1]; c2=cadj[i][2];
    p1=cp[c1]; p2=cp[c2];
    if ( j>i || j1!=0) {
      if (rescaled && !leftfil && PERX(j1)<0) goto SKIPCARC;
      if (rescaled && !botfil && PERY(j1)<0) goto SKIPCARC;
      trans(vx[j],vy[j],j1,&x2,&y2);
      alpha=carcangle(x1,y1,x2,y2,p1,p2);
      straight= !arccentre(x1,y1,x2,y2,alpha,&xc,&yc);
      if (straight) {
        if (rescaled) {
          x1=xorigin+scale*x1; y1=yorigin+scale*y1;
          x2=xorigin+scale*x2; y2=yorigin+scale*y2;
        }
        mglline(x1,y1,x2,y2);
      }
      else {
        r=linlen(x1,y1,xc,yc);
        th=linangle(xc,yc,x1,y1);
        if (rescaled) {
          r *= scale;
          xc=xorigin+scale*xc; yc=yorigin+scale*yc;
        }
        th1=th-alpha;
        if (th1>th) {
          temp=th; th=th1; th1=temp;
        }
        mglarc(xc,yc,r,th1,th);
      }
    }
SKIPCARC:
    /* Case of border arcs */
    x1=vx[i]; y1=vy[i];
    for (k=1; k<3; k++) {
      j=vnbr[i][k]; j1=(vper[i][k] & PERMASK);
      if ( j>i || j1!=0) {
        if (rescaled && !leftfil && PERX(j1)<0) goto SKIPBARC;
        if (rescaled && !botfil && PERY(j1)<0) goto SKIPBARC;
        trans(vx[j],vy[j],j1,&x2,&y2);
        la=larc(i,k); b=cadj[i][0];
        if (k==1) alpha=barcangle(x1,y1,x2,y2,p2,bp[b],la,&arcbrk);
        if (k==2) alpha=barcangle(x1,y1,x2,y2,bp[b],p1,la,&arcbrk);
        if (arcbrk) {
          plerror("unplottable border arc in mglfoam");
          goto SKIPBARC;
        }
        straight= !arccentre(x1,y1,x2,y2,alpha,&xc,&yc);
        if (straight) {
          if (rescaled) {
            x1=xorigin+scale*x1; y1=yorigin+scale*y1;
            x2=xorigin+scale*x2; y2=yorigin+scale*y2;
          }
          mglline(x1,y1,x2,y2);
        }
        else {
          r=linlen(x1,y1,xc,yc);
          th=linangle(xc,yc,x1,y1);
          if (rescaled) {
            r *= scale;
            xc=xorigin+scale*xc; yc=yorigin+scale*yc;
          }
          th1=th-alpha;
          if ((twopimod(th-th1)>PI && !la) || (twopimod(th-th1)<PI && la)) {
            temp=th; th=th1; th1=temp;
          }
          mglarc(xc,yc,r,th1,th);
        }
      }
SKIPBARC: ;
    }
  }
  /* Now draw the isolated bubbles, if any */
  for (ii=0; ii<nbub; ii++) {
    i=bublist[ii];
    xc=cx[i]; yc=cy[i]; bubinbox(&xc,&yc);
    r=BRADIUS(cp[i],bpav);
    if (rescaled) { xc=xorigin+scale*xc; yc=yorigin+scale*yc; r *= scale; }
    mglarc(xc,yc,r,0.0,2.0*PI);
  }
  if (hatch) {
    for (i=0; i<onb; i++) found[i]=FALSE;
    for (ii=0; ii<nv; ii++) {
      i=vlist[ii];
      if (!found[b=cadj[i][0]]) {
        found[b]= TRUE;
        mglhatch(i,theta,xorigin,yorigin,scale);
      }
    }
  }  
	if (centroid) {
		for(i=0;i<nc;i++){
			mglpoint(cxcent[i],cycent[i]);
		}
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *    M G L B O R D ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	mglbord(short i)
 *
 *	Arguments:	i	= index of vertex lying on a Plateau border
 *
 *	Return value:	none
 *
 *	Action:		Plot the Plateau border adjacent to vertex `i'.
 *			This can be used to plot interesting large Plateau
 *			borders -- in particular as percolation is approached.
 *
 *****************************************************************************/
void mglbord(i)
short i;
{
  short j, j1;
  short perconcat();
  REAL x1, y1, x2, y2, p1, pb, alpha, xc, yc, r, th, th1, temp;
  REAL barcangle(), carcangle(), linangle(), linlen(), twopimod();
  boolean straight, arcbrk, la, arccentre(), larc();
  void mglline(), mglarc(), trans(), bubinbox();
  void plerror();
  x1=vx[i]; y1=vy[i]; pb=bp[cadj[i][0]];
  j=i; j1=0;
  do {
    p1=cp[cadj[j][2]]; la=larc(j,1);
    j1=perconcat(j1,vper[j][1]); j=vnbr[j][1];
    trans(vx[j],vy[j],j1,&x2,&y2);
    alpha=barcangle(x1,y1,x2,y2,p1,pb,la,&arcbrk);
    if (arcbrk) {
      plerror("unplottable border arc in routine mglbord"); continue;
    }
    straight= !arccentre(x1,y1,x2,y2,alpha,&xc,&yc);
    if (straight) {
      mglline(x1,y1,x2,y2);
    }
    else {
      r=linlen(x1,y1,xc,yc);
      th=linangle(xc,yc,x1,y1);
      th1=th-alpha;
      if ((twopimod(th-th1)>PI && !la) || (twopimod(th-th1)<PI && la)) {
        temp=th; th=th1; th1=temp;
      }
      mglarc(xc,yc,r,th1,th);
    }
    x1=x2; y1=y2;
  } while (j != i);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *  M G L H A T C H ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	mglhatch(short i, REAL theta, REAL xorigin,
 *				 REAL yorigin, REAL scale)
 *
 *	Arguments:	i	= index of vertex which lies on the Plateau
 *				border
 *			theta	= sloping angle of hatching lines
 *			(xorigin, yorigin) = centre point of network
 *				(ignored if `scale == 0')
 *			scale	= scaling factor
 *				(defaults to 1.0 if scale == 0)
 *
 *	Return value:	none
 *
 *	Action:		Fills the Plateau border (specified by giving a vertex
 *			`i' which lies on it) with hatching lines.  The slope
 *			of the hatching lines is given by `theta'.  If the
 *			centre of the network has been shifted then it is
 *			given by (xorigin, yorigin) and all points are 
 *			rescaled by the factor `scale'.
 *
 *****************************************************************************/
void mglhatch(i,theta,xorigin,yorigin,scale)
short i;
REAL theta, xorigin, yorigin, scale;
{
  REAL ptx[MSEGMENTS], pty[MSEGMENTS];
  REAL xx, yy, d, dc, r, dmin, dmax, da, db, xa, ya, xb, yb,
       xp, yp, xq, yq, x1, y1, x2, y2, h, alpha, p1, pb,
       xc1, yc1, th, thp, thq, wid2, hgt2, t;
  short i1, ii, ii1, j, j1, je, l, nwh, wheel[MWHEEL][2], npt, b, n;
  boolean rescaled, la, arcbrk, straight, arccentre(), larc();
  REAL barcangle(), linangle(), twopimod(), intersect_cmp(), fsign();
  short perconcat();
  void setlarc(), trans(), intersect_ins(), mglline();
  void plerror();
  rescaled = (scale!=0);
  /* Define a wheel... */
  b=cadj[i][0];  pb=bp[b];
  nwh=0;
  wheel[0][0]=i;
  i1=vper[i][1];  wheel[0][1]=(i1 & LARC);
  nwh++;
  je=i;
  while ((j=vnbr[i][1]) != je) {
    j1=perconcat(i1,vper[j][1]);
    setlarc(&j1,larc(j,1));
    wheel[nwh][0]=j;
    wheel[nwh][1]=i1;
    nwh++;
    i=j;
    i1=j1;
  }
  /* now find out where the first hatched line goes */
  xx=cos(PI*theta/180.0);
  yy=sin(PI*theta/180.0);
  dmin=1.0e6;  dmax= -1.0e6;
  for (l=0; l<nwh; l++) {
    i=wheel[l][0];  i1=wheel[l][1];
    trans(vx[i],vy[i],i1,&x1,&y1);
    d= -x1*yy + y1*xx;
    dmin=min(dmin,d);
    dmax=max(dmax,d);
  }
  /* now loop over all hatched lines */
  h= HATCHSPACING;
  for (d=dmin+h; d<dmax; d += h) {
    /* create the list of points where the 'hatch-line' */
    /* intersects the Plateau border                    */
    npt=0;
    for (l=0; l<nwh; l++) {
      i=wheel[l][0];
      i1=wheel[l][1];
      ii=wheel[(l+1)%nwh][0];
      ii1=wheel[(l+1)%nwh][1];
      trans(vx[i],vy[i],i1,&x1,&y1);
      trans(vx[ii],vy[ii],ii1,&x2,&y2);
      la=larc(i,1);
      p1=cp[cadj[i][2]];
      alpha=barcangle(x1,y1,x2,y2,p1,pb,la,&arcbrk);
      straight= !arccentre(x1,y1,x2,y2,alpha,&xc1,&yc1);
      if (!straight) {
        dc= -xc1*yy + yc1*xx;
        r=BRADIUS(p1,pb);
        if (dc-r<d && d<dc+r) {
          /* intersection with the circle */
          da=d-dc;
          xa= -da*yy;  ya=da*xx;
          db=sqrt(r*r-da*da);
          xb=db*xx;  yb=db*yy;
          xp=xc1+xa+xb;  yp=yc1+ya+yb;
          xq=xc1+xa-xb;  yq=yc1+ya-yb;
          th=linangle(xc1,yc1,x2,y2);
          thp=linangle(xc1,yc1,xp,yp);
          thq=linangle(xc1,yc1,xq,yq);
          if (twopimod(thp-th)<alpha) {
            /* a point of intersection is found */
            n=0;
            while (n<npt && intersect_cmp(xp,yp,ptx[n],pty[n])<0.0) n++;
            intersect_ins(&npt,ptx,pty,n,xp,yp);
          }
          if (twopimod(thq-th)<alpha) {
            /* a point of intersection is found */
            n=0;
            while (n<npt && intersect_cmp(xq,yq,ptx[n],pty[n])<0.0) n++;
            intersect_ins(&npt,ptx,pty,n,xq,yq);
          }
        }
      }
    }
    /* impose vertical boundaries... */
    if (npt%2!=0) plerror("mglhatch: oddity error!");
    for (n=0; n<npt; n += 2) {
      wid2=0.5*boxwid;
      xp=ptx[n];  yp=pty[n];
      xq=ptx[n+1];  yq=pty[n+1];
      if (fabs(xp)>wid2 && fabs(xq)>wid2) {
        t=fsign(xp)*boxwid;
        ptx[n] -= t;  ptx[n+1] -= t;
      }
      if (fabs(xp)>wid2 && fabs(xq)<=wid2) {
        t=fsign(xp);
        xa=t*wid2;
        ya=yp+(yq-yp)*(xa-xp)/(xq-xp);
        intersect_ins(&npt,ptx,pty,n+1,xa,ya);
        intersect_ins(&npt,ptx,pty,n+1,xa,ya);
        t *= boxwid;
        ptx[n] -= t;  ptx[n+1] -= t;
        n += 2;
      }
      if (fabs(xp)<=wid2 && fabs(xq)>wid2) {
        t=fsign(xp);
        xa=t*wid2;
        ya=yp+(yq-yp)*(xa-xp)/(xq-xp);
        intersect_ins(&npt,ptx,pty,n+1,xa,ya);
        intersect_ins(&npt,ptx,pty,n+1,xa,ya);
        t *= boxwid;
        ptx[n+2] -= t;  ptx[n+3] -= t;
        n += 2;
      }
    }
    /* impose horizontal boundaries... */
    for (n=0; n<npt; n += 2) {
      hgt2=0.5*boxhgt;
      xp=ptx[n];  yp=pty[n];
      xq=ptx[n+1];  yq=pty[n+1];
      if (fabs(yp)>hgt2 && fabs(yq)>hgt2) {
        t=fsign(yp)*boxhgt;
        pty[n] -= t;  pty[n+1] -= t;
      }
      if (fabs(yp)>hgt2 && fabs(yq)<=hgt2) {
        t=fsign(yp);
        ya=t*hgt2;
        xa=xp+(xq-xp)*(ya-yp)/(yq-yp);
        intersect_ins(&npt,ptx,pty,n+1,xa,ya);
        intersect_ins(&npt,ptx,pty,n+1,xa,ya);
        t *= boxhgt;
        pty[n] -= t;  pty[n+1] -= t;
        n += 2;
      }
      if (fabs(yp)<=hgt2 && fabs(yq)>hgt2) {
        t=fsign(yp);
        ya=t*hgt2;
        xa=xp+(xq-xp)*(ya-yp)/(yq-yp);
        intersect_ins(&npt,ptx,pty,n+1,xa,ya);
        intersect_ins(&npt,ptx,pty,n+1,xa,ya);
        t *= boxhgt;
        pty[n+2] -= t;  pty[n+3] -= t;
        n += 2;
      }
    }
    /* now plot the list of segments */
    for (n=0; n<npt; n += 2) {
      xp=ptx[n];  yp=pty[n];
      xq=ptx[n+1];  yq=pty[n+1];
      if (rescaled) {
        xp=xorigin+scale*xp;
        yp=yorigin+scale*yp;
        xq=xorigin+scale*xq;
        yq=yorigin+scale*yq;
      }
      mglline(xp,yp,xq,yq);
    }
  }
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *   I N T E R S E C T _ I N S ( )  * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	intersect_ins(short *npt, REAL ptx[], REAL pty[],
 *				      short n, REAL x1, REAL y1)
 *
 *	Arguments:	*npt	= size of following arrays
 *			(ptx[], pty[]) = coordinates of line segments
 *			n	= index where new point is to be inserted
 *			(x1, y1) = coordinates of point to be inserted
 *
 *	Return value:	none
 *
 *	Action:		Inserts the point (x1, y1) at index n in the arrays
 *			ptx[] and pty[] (pushing up the other elements).
 *			Then increments the total number of points *npt.
 *
 *****************************************************************************/
void intersect_ins(npt,ptx,pty,n,x1,y1)
short *npt, n;
REAL ptx[], pty[], x1, y1;
{
  short k;
  for (k= *npt; k>n; k--) {
    ptx[k]=ptx[k-1];  pty[k]=pty[k-1];
  }
  ptx[n]=x1;  pty[n]=y1;
  *npt += 1;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *   I N T E R S E C T _ C M P ( )  * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	intersect_cmp(REAL x1, REAL y1, REAL x2, REAL y2)
 *
 *	Arguments:	(x1, y1) = p1 -- first 2-d point
 *			(x2, y2) = p2 -- second 2-d point
 *
 *	Return value:	Positive if `p1 > p2'; else Negative.
 *
 *	Action:		Defines a convenient ordering of points in the
 *			2-dimensional plane.
 *
 *****************************************************************************/
REAL intersect_cmp(x1,y1,x2,y2)
REAL x1, y1, x2, y2;
{
  void plerror();
  if (x1!=x2) { return(x1-x2); }
  if (y1!=y2) { return(y1-y2); }
  plerror("intersect_cmp: you should not get here!");
  return(0.0);
}
