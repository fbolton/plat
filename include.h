#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "machine.h"

#ifdef IRIS
#include <gl.h>
#endif

#ifdef X_WINDOWS
#include <Xm/Xm.h>
#endif

#ifndef boolean
#define boolean int
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef max
#define max(a,b) (((a)>(b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b) (((a)<(b)) ? (a) : (b))
#endif

/* The following three parameters MVERT, MBORD and MCELL all bear exact
ratios to each other (as shown) and *must* be changed together. */

#define MVERT 30000  /* Max no. of vertices:  MVERT=6*MCELL */

#define MBORD 10000  /* Max no. of Plateau borders:  MBORD=2*MCELL */

#define MCELL 5000  /* Max no. of cells:  Fixes the above parameters. */

#define MSTRING 256  /* Max length of a command string */

#define MINFO 16  /* Update this whenever you add another type of
information to the list in 'plat.c' */

#define VORRAD 0.2  /* The radius used for the initial
Voronoi triangulation search circle.  Roughly in units of a typical
cell diameter */

#define HARDRFRAC 0.65  /* This specifies the minimum distance between any
two points in the voronoi network.  It is specified as a fraction of the 
hexagonal close packing value, thus '1.0' is the max possible! */

#define VORSEED 1  /* Can be any +ve int (a random number seed) */

#define MINBFRAC 0.05  /* Min value of the 'radius' of a Plateau
border which is decorated onto a vertex during initialisation.
It is expressed as a fraction of the typical radius of a cell. */

#define MAXBFRAC 0.05  /* Max value of the 'radius' ... (ditto) */

#define PHIFRACTION 0.2  /* Gives the fractional change in size of a
typical PB per step of the routine 'makefraction', which sets the
area-fraction of the foam. */

#define VDAMP 0.9  /* Damping factor used to relax an ordinary vertex */

#define VVDAMP 0.5 /* Damping factor used in relaxing a vertex which is
adjacent to a PB with more than three sides */

#define PDAMP 0.9  /* Damping factor used for pressure equilibration */

#define BPDAMP 0.5  /* Damping factor used for equilibration of PB pressure */

#define PHIDAMP 0.5  /* Damping for the routine 'phirelax' used to alter the
volume fraction of the foam */

#define BPRELAX 1.0  /* This determines the rate at which a border pressure
'bp[i]' relaxes back to the averate value 'bpav' (in 'normalizep()') */

#define PI 3.1415926536

#define DFLTLENGTHDELTA 0.0005  /* 'delta' used when doing numerical
differentiation wrt distance. Expressed as a fraction of the typical
PB radius. */

#define DFLTPRESSUREDELTA 0.0001  /* 'delta' used when doing numerical
differentiation wrt pressure. The value given here results in a change
of about 1 in 10e4 to the cell pressure. It is expressed as a fraction of the
pressure of an isolated bubble with the average area of the network.
This figure is also used as a
guide for when a pressure increment is considered so small as to be
negligeable */

#define DFLTDIFFUSERATE 0.01  /* this gives the typical percentage
of a cell's area which is diffused */

#define DFLTAREASUP 0.05  /* this gives a typical value for the allowed
error in the area of a cell. */

#define MINAREASUP 0.001 /* the min fractional error in the area of a cell */

#define DFLTBAREASUP 1.0 /* this gives the allowed error in the area of
a Plateau border */

#define DFLTMAXDV 0.05  /* sets the largest average increment given to a
triangular Plateau border - expressed as a fraction of the
typical cell radius */

#define DFLTMAXDVV 0.01 /* sets the largest allowed distortion of a Plateau
border - given as a fraction of the typical border radius */

#define MAXITER 1000  /* sets the max. no. of iterations per equilibration */

#define DFLTEQUILSUP 0.02  /* taken as the usual convergence value
for the 'sup' which is returned by 'equil()'. It is expressed as a
fraction of the PB size (NB See DFLTMINVVLEN) */

#define DFLTMINVVLEN 0.01  /* Min vertex to vertex length used to determine
seperation/coalescence - expressed as a fraction of the PB radius. Must
be of the same order as DFLTEQUILSUP. */

#define DFLTFILMWID 0.0  /* Used by 'bpinch()' to determine if a film is
thin enough to undergo a spontaneous pinch */

#define COSMINANG 0.999999  /* This determines the smallest angle deviation
from 'pi' before angle 'popping' occurs */

/* Sets the max no. of failed indices recorded in Equil */
#define MFAILEDINDEX 20

#define MWHEEL 30 /* ~max no. of sides to a cell, ~max no. of sides
on a Plateau border */

#define MPRECISION 1.0e-4  /* precision of calculations (sort of!) */

#define MELOST 30  /* max no. of edge losses which may be recorded
per network relaxation */

#define STUCK 2  /* used as part of a flag in the pressure equilibration loop*/

#define TOFILE 2 /* used to indicate that info will be saved on file */

#define TOSCREEN 1 /* info to be send to 'stdout' */

#define MSEGMENTS 1000  /* gives the max. no. of points in a hatch-line */

#define SET_HATCHING 0

#define HATCHSPACING 0.001  /* given as a fraction of box size */

/* Include these lines to use an IRIS (Silicon Graphics) terminal */
#ifdef IRIS
#define IRISXO 300.0 /*700.0*/ /* usually 300.0 */
#define IRISYO 300.0 /*-400.0*/
#define IRISXS 450.0 /*3000.0*/ /* usually 450.0 */
#define IRISYS 450.0 /*3000.0*/
#define IRISPREFLX 400
#define IRISPREFRX 1000
#define IRISPREFBY 200
#define IRISPREFTY 800
#endif

/* Include these lines to use the RD-GL Plotter */
/* ...or the O'Reilly HP-GL Plotter             */
/* (Hewlett-Packard graphics language) */
#if defined(HPGL) || defined(ORIHPGL)
#define HPGLXO 8000.0
#define HPGLYO 5500.0
#define HPGLXS 6000.0
#define HPGLYS 6000.0
#endif

/* Include these lines to generate a Tektronix file */
#ifdef TEK
#define TEKXO 365.0
#define TEKYO 365.0
#define TEKXS 400.0
#define TEKYS 400.0
#define TEKARCNUM 10   /* Number of segments to an arc. */
#endif

/* Include these lines to generate Postscript graphics */
#ifdef POSTSCRIPT
#define POSTXO 1200  /* usually 300 */
#define POSTYO 1600  /* usually 400 */
#define POSTXS 1800  /* usually 450 */
#define POSTYS 1800  /* usually 450 */
#define POSTARCNUM  10
#define POSTLINEWID 1.0  /* measured in 'mm' */
#endif

/* Include these lines to use X11 Release 5 and Motif1.2 graphics */
#ifdef X_WINDOWS
#define X_WIN_XO	300
#define X_WIN_YO	300
#define X_WIN_XS	600
#define X_WIN_YS	600
#endif

/* These macros are physically very important.  They define the */
/* relation between arc radii and pressure differences. */
#define BRADIUS(p1,p2) ((1.0)/((p1)-(p2)))
#define CRADIUS(p1,p2) ((2.0)/((p1)-(p2)))

/* These macros, in large part, implement the boundary conditions */
#define PERFN(ix,iy) perfn((ix),(iy))
#define PERX(j) perx(j)
#define PERY(j) pery(j)
#define PERINCX(j) perincx(j)
#define PERDECX(j) perdecx(j)
#define PERINCY(j) perincy(j)
#define PERDECY(j) perdecy(j)

#define PERMASK 0x0ff /* masks that part of an index which gives the periodic
index */

#define LARC 0x100 /* used to set the 'large arc' flag bit in
a periodic index */

/* Global Function Declarations */
short perfn(), perx(), pery(), perincx(), perdecx(), perincy(), perdecy();
void plerror();

/* Voronoi variables */
extern boolean waspivot[MCELL];
extern short triang[MBORD][3][2], vorvnbr[MBORD][3], vorcadj[MBORD][3],
      vorvper[MBORD][3];
extern REAL cx[MCELL], cy[MCELL];

/* Variables to facilitate X windows */
#ifdef X_WINDOWS
struct bitmap_struct {
  GC draw_gc, undraw_gc;
  Pixmap froth_bitmap;
  Widget frothPicture;
  String filename;
  Dimension pixmap_width, pixmap_height;
};
extern struct bitmap_struct bitmap;
#endif

/* Variables to control graphics */
extern boolean mgl_hpgl_flag, mgl_tek_flag, mgl_ps_flag;
extern char mgl_hpgl_filename[256], mgl_tek_filename[256], mgl_ps_filename[256];


/* General variables */
extern short knbr[3];
extern REAL vx[MVERT], vy[MVERT], cp[MCELL], carea[MCELL], darea[MCELL],
       dvx[MVERT], dvy[MVERT], dcp[MCELL];
extern REAL bpav, bp[MBORD], barea[MBORD];
extern short vlist[MVERT], clist[MCELL], blist[MBORD], bublist[MCELL],
             vnbr[MVERT][3], vper[MVERT][3], cadj[MVERT][3];
extern short nbsides[MBORD], ncsides[MCELL], iel[MELOST];
extern short nv, nc, nb, nbub, onv, onc, onb, nel;
extern int elosscount, bpinchcount;
extern REAL tfoam, boxwid, boxhgt, netenergy, henckyeps, volfrac, minenergy;
extern REAL lengthdelta, pressuredelta, equilsup, filmwid, minvvlen,
            diffuserate, bprelax, areasup, bareasup, minbfrac, maxbfrac,
            maxdv, maxdvv, cosminang, bscale, cscale;
extern REAL svdamp, svvdamp, spdamp, sbpdamp, sbprelax, sminbfrac, smaxbfrac,
            slengthdelta, spressuredelta, sdiffuserate, sareasup,
            sbareasup, smaxdv, smaxdvv, sequilsup, sminvvlen, sfilmwid,
            scosminang, srfrac, snotopol, sminareasup, smaxiter, svorseed;
extern boolean notopol;
extern char info_tok[][20];
extern boolean info_list[MINFO], info_wlist[MINFO], foamlike, ginteract;
