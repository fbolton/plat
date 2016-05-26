/*
 * This file implements an X-windows driven interface to the 'plat' program.
 * It is based on OSF/Motif version 1.2 and on the X11 release 5.  The program
 * makes calls to 'motif', 'X intrinsics' and 'Xlib' all of which are
 * localised in this file.
 *
 * To enable this X interface just set the 'X_WINDOWS' flag in the file
 * 'machine.h'.
 *
 * F. Bolton 17 Feb 1995.
 *
 */

/* Standard Motif include files and other includes. */
#include "include.h"
#include <string.h>
#include <ctype.h>
#include <Xm/RepType.h>

/*
 * Public header files for Widgets used in this file.
 */
#include <Xm/MainW.h>
#include <Xm/RowColumn.h>
#include <Xm/PanedW.h>
#include <Xm/Form.h>
#include <Xm/Frame.h>
#include <Xm/PushB.h>
#include <Xm/CascadeB.h>
#include <Xm/MessageB.h>
#include <Xm/Command.h>
#include <Xm/DrawingA.h>

/*
 * Defined macros: will ultimately be moved to 'include.h'
 */

#define X_PIXMAP_WIDTH		594
#define X_PIXMAP_HEIGHT		594

/* 
 * X-WINDOW GLOBAL VARIABLES
 */

XtAppContext app_context;
struct bitmap_struct bitmap;
static String command_string;
static char   command_previous[256];
static boolean x_window_mode, xuser_show_equil, xuser_pause, xuser_stop;

/*
 * Callback functions defined here.
 */

/*ARGSUSED*/
void Quit(Widget w, XtPointer client_data, XtPointer call_data)
{
  exit(0);
}

/*ARGSUSED*/
void UpdateScreen(Widget w, XtPointer client_data, XtPointer call_data)
{
  if (DefaultDepthOfScreen(XtScreen(w)) == 1)
    XCopyArea(XtDisplay(w),
	bitmap.froth_bitmap,
	XtWindow(w),
	DefaultGCOfScreen(XtScreen(w)), 0, 0,
	bitmap.pixmap_width,
	bitmap.pixmap_height,
	0, 0);
  else
    XCopyPlane(XtDisplay(w),
	bitmap.froth_bitmap,
	XtWindow(w),
	DefaultGCOfScreen(XtScreen(w)), 0, 0,
	bitmap.pixmap_width,
	bitmap.pixmap_height,
	0, 0, 1);
}

/*ARGSUSED*/
void ResetPlot(Widget w, XtPointer client_data, XtPointer call_data)
{
  void mglreset(), foamplot();
  Widget wtop;

  mglreset();
  foamplot(0);
}

/*ARGSUSED*/
void CommandUpdate(Widget w, XtPointer client_data, XtPointer call_data)
{
  XmString compound;

  XtVaGetValues(w, XmNcommand, &compound, NULL);
  (void) XmStringGetLtoR(compound, XmFONTLIST_DEFAULT_TAG,
 	 &command_string);
  XmStringFree(compound);
  if (strlen(command_string) > 0)
    strcpy(command_previous,(char *) command_string);
}

/*ARGSUSED*/
void PlatCommand(Widget w, XtPointer client_data, XtPointer call_data)
{
  void plerror(), setconstants(), trans(), vorspray();
  void mglinit(), mglline(), mglclose(), mglreset();
  void vfoamplot(), foamplot(), polyplot(), voronoi(), hexnet(),
       vortrans(), ptbox(), hencky(), s_execute();
  short i, j, k, nestlevel, slen;
  REAL x1, y1, x2, y2, d, eps;
  float eps0, f;
  char s[256], c;
  boolean sgood;

  strcpy(s, command_previous);
  nestlevel=0; slen=0; sgood=TRUE;
  xuser_show_equil = -1;  /* Neutral setting for equilibration plotting */
  xuser_stop = FALSE;
  xuser_pause = FALSE;
  do {
    if (sgood) {
      slen=strlen(s);
      /* First check if the execution of the string ought to be suspended */
      /* pending further input */
      nestlevel=0;
      for (i=0; i<slen; i++) {
        switch(s[i]) {
        case '{' : nestlevel++; break;
        case '}' : nestlevel--;
                   if (nestlevel<0) {
                     plerror("too many right brackets");
                     sgood=FALSE; i=slen;
                   }
                   break;
        case '\\' :
          s[i]=' ';
          if (i!=slen) {
            plerror("illegal embedded backslash character \'\\\'");
            sgood=FALSE; i=slen;
          }
          break;
        }
      }
    }
    s[slen]='\n'; slen++; s[slen]='\0';
  } while (nestlevel>0 && sgood);
  if (sgood) s_execute(s,0);
  fflush(stdout);
}

/*ARGSUSED*/
static void
ShowEquil(Widget w, XtPointer client_data, XtPointer call_data)
{
  if (xuser_show_equil== -1 || xuser_show_equil==0)
    xuser_show_equil = 1;	/* Emphatically *do* plotting */
  else if (xuser_show_equil==1)
    xuser_show_equil = 0;	/* Emphatically *don't* do plotting */
}

/*ARGSUSED*/
static void
Pause(Widget w, XtPointer client_data, XtPointer call_data)
{
  if (xuser_pause) 	xuser_pause = FALSE;
  else 			xuser_pause = TRUE;
}

/*ARGSUSED*/
static void
Stop(Widget w, XtPointer client_data, XtPointer call_data)
{
  xuser_stop = TRUE;
}

/*
 * Action Functions Defined Here
 */

/*ARGSUSED*/
static void
ZoomIn(Widget w, XButtonEvent *event, String *params, Cardinal *num_params)
{
  void zoom(), foamplot();

  zoom(1, event->x, event->y);
  foamplot(0);
}

/*ARGSUSED*/
static void
ZoomCentre(Widget w, XButtonEvent *event, String *params, Cardinal *num_params)
{
  void zoom(), foamplot();

  zoom(2, event->x, event->y);
  foamplot(0);
}

/*ARGSUSED*/
static void
ZoomOut(Widget w, XButtonEvent *event, String *params, Cardinal *num_params)
{
  void zoom(), foamplot();

  zoom(3, event->x, event->y);
  foamplot(0);
}

/*
 * GLOBAL 'PLAT' VARIABLES AND VARIABLES FOR COMMAND DRIVER
 */

/* Any changes made to these Global Variables must be reflected in the file
'include.h' to avoid a general state of chaos. */

/* Voronoi variables */
boolean waspivot[MCELL];
short triang[MBORD][3][2], vorvnbr[MBORD][3], vorcadj[MBORD][3],
      vorvper[MBORD][3];
REAL cx[MCELL], cy[MCELL];

/* Variables to control graphics */
boolean mgl_hpgl_flag= FALSE, mgl_tek_flag= FALSE, mgl_ps_flag= FALSE;
char mgl_hpgl_filename[256], mgl_tek_filename[256], mgl_ps_filename[256];

/* General variables */
short knbr[3] = {0, 2, 1};
REAL vx[MVERT], vy[MVERT], cp[MCELL], carea[MCELL], darea[MCELL], dvx[MVERT],
      dvy[MVERT], dcp[MCELL], cxcent[MCELL], cycent[MCELL];
REAL bpav=0.0, bp[MBORD], barea[MBORD];
short vlist[MVERT], clist[MCELL], blist[MBORD], bublist[MCELL], vnbr[MVERT][3],
      vper[MVERT][3], cadj[MVERT][3];
short nbsides[MBORD], ncsides[MCELL], iel[MELOST];
short nv, nc, nb, nbub=0, onv, onc, onb, nel;
int elosscount, bpinchcount;
REAL tfoam=0.0, boxwid=1.0, boxhgt=1.0, netenergy, henckyeps, volfrac=0.0, minenergy;
REAL lengthdelta, pressuredelta, equilsup, filmwid, minvvlen,
      diffuserate, bprelax, areasup, bareasup, minbfrac, maxbfrac, maxdv,
      maxdvv, cosminang, bscale, cscale;
REAL svdamp= VDAMP, svvdamp= VVDAMP, spdamp= PDAMP, sbpdamp= BPDAMP,
      sbprelax= BPRELAX, sminbfrac= MINBFRAC, smaxbfrac= MAXBFRAC,
      slengthdelta= DFLTLENGTHDELTA, spressuredelta= DFLTPRESSUREDELTA,
      sdiffuserate= DFLTDIFFUSERATE, sareasup= DFLTAREASUP,
      sbareasup= DFLTBAREASUP, smaxdv= DFLTMAXDV, smaxdvv= DFLTMAXDVV,
      sequilsup= DFLTEQUILSUP, sminvvlen= DFLTMINVVLEN,
      sfilmwid= DFLTFILMWID, scosminang= COSMINANG, srfrac= HARDRFRAC,
      snotopol= 0.0, sminareasup= MINAREASUP, smaxiter= MAXITER,
      svorseed= (REAL) VORSEED;
/*                      * Infolist */
/* NB: Whenever this list is updated, MINFO must be changed in 'include.h' */
char info_tok[][20]={"muside", "rhoside", "muarea", "arean", "bside",
     "netenergy", "henckyeps", "elosscount", "bpinchcount", "phi",
     "nbar", "edgefrac", "z", "arootbar", "minenergy", "adjmat", "centroids"};
boolean info_list[MINFO], info_wlist[MINFO], foamlike, ginteract=FALSE,
        notopol=FALSE;

main(int argc, char **argv)
{
  /*
   * 'plat' Initialisation Routines:
   */
  void mglreset();
  void plat_command_driver();

  /*
   * 'X' Window Variables...
   */
  void SetUpPixmap(Widget);
  void Test();
  Widget topLevel, mainWindow, menuBar, frame;
  Widget fileButton, fileMenu, quit, pauseButton, stopButton,
	 resetPlotButton, showEquilButton, platCommand, workSpace;
  Widget leftButtonCol, rightButtonCol;
  Widget rowForm;
  static char transTable[] = 
    "<Btn1Down>:	ZoomIn()\n\
     <Btn2Down>:	ZoomCentre()\n\
     <Btn3Down>:	ZoomOut()";
  static XtActionsRec actions[] = {
    {(String) "ZoomIn", (XtActionProc) ZoomIn},
    {(String) "ZoomOut", (XtActionProc) ZoomOut},
    {(String) "ZoomCentre", (XtActionProc) ZoomCentre}
  };


  /*
   * Local Variables
   */
  int i;
  
  /*
   * Perform 'plat' initialisation
   */
  mglreset();
  nc=MCELL;
  for (i=0; i<MINFO; i++) {
    info_list[i]=FALSE; info_wlist[i]=FALSE;
  }
  x_window_mode = FALSE;	/* initially 'X' windows not used */
  xuser_show_equil = -1;
  xuser_pause  = FALSE;
  xuser_stop   = FALSE;

  /* 
   * Call the Command Driver
   */
  plat_command_driver();

  /*
   * Enter the 'X' window interface at this point
   * (occurs if the user selects option 'g' from the command driver)
   */
  XtSetLanguageProc(NULL, (XtLanguageProc)NULL, NULL);

  topLevel = XtVaAppInitialize(
	  &app_context,		/* Application context */
	  "XPlat",		/* application class name */
	  NULL, 0, 		/* command line option list */
	  &argc, argv,		/* command line args */
	  NULL,			/* for missing app-defaults file */
	  NULL);		/* terminate varargs list */
  
  /* create the main window */
  mainWindow = XtVaCreateManagedWidget(
	  "mainWindow",		/* widget name */
	  xmMainWindowWidgetClass,
	  topLevel,		/* parent widget */
	  NULL);

  /* register converter for setting tearoff menus from resource files */
  XmRepTypeInstallTearOffModelConverter();

  /* create menubar along top inside of main window */
  menuBar = XmCreateMenuBar(
	  mainWindow,		/* parent widget */
	  "menuBar",		/* widget name */
	  NULL, 0);		/* no arguments needed */
  XtManageChild(menuBar);

  workSpace = XtVaCreateManagedWidget(
	  "workSpace",
	  xmFormWidgetClass,
	  mainWindow,
	  NULL);

  /*
   * SET UP THE BITMAP AND DRAWING AREA WIDGET
   */

  bitmap.pixmap_height = X_PIXMAP_HEIGHT;
  bitmap.pixmap_width  = X_PIXMAP_WIDTH;
  /* Create pixmap, GC and draw initial contents into pixmap */
  SetUpPixmap(topLevel);

  frame = XtVaCreateManagedWidget(
	  "frame",
	  xmFrameWidgetClass,
	  workSpace,
	  NULL);

  bitmap.frothPicture = XtVaCreateManagedWidget(
	  "frothPicture",
	  xmDrawingAreaWidgetClass,
	  frame,
	  XmNtranslations, XtParseTranslationTable(transTable),
	  XmNwidth, bitmap.pixmap_width,
	  XmNheight, bitmap.pixmap_height,
	  NULL);

  XtAddCallback(bitmap.frothPicture, XmNexposeCallback, UpdateScreen, 0);

  /*
   * SET UP AN ARRAY OF FOUR LARGE BUTTONS IN FRAMES
   */

  rowForm = XtVaCreateManagedWidget(
	  "rowForm",
	  xmFormWidgetClass,
	  workSpace,
	  NULL);

  leftButtonCol = XtVaCreateManagedWidget(
	  "leftButtonCol",
	  xmRowColumnWidgetClass,
	  rowForm,
	  NULL);

  rightButtonCol = XtVaCreateManagedWidget(
	  "rightButtonCol",
	  xmRowColumnWidgetClass,
	  rowForm,
	  NULL);

  resetPlotButton = XtVaCreateWidget(
	  "resetPlotButton",
	  xmPushButtonWidgetClass,
	  leftButtonCol,
	  NULL);

  showEquilButton = XtVaCreateWidget(
	  "showEquilButton",
	  xmPushButtonWidgetClass,
	  leftButtonCol,
	  NULL);

  pauseButton = XtVaCreateWidget(
	  "pauseButton",
	  xmPushButtonWidgetClass,
	  rightButtonCol,
	  NULL);

  stopButton = XtVaCreateWidget(
	  "stopButton",
	  xmPushButtonWidgetClass,
	  rightButtonCol,
	  NULL);

  XtManageChild(resetPlotButton);
  XtManageChild(showEquilButton);
  XtManageChild(pauseButton);
  XtManageChild(stopButton);

  /*
   * MAKE THE COMMAND WIDGET...
   */

  platCommand = XtVaCreateManagedWidget(
	  "platCommand",
          xmCommandWidgetClass,
          mainWindow,
          NULL);

  /*
   * CREATE FILE MENU AND CHILDREN
   */

  /* create the file button in the menubar */
  fileButton = XtVaCreateManagedWidget(
	  "fileButton",
	  xmCascadeButtonWidgetClass,
	  menuBar,
	  NULL);

  /* create menu */
  fileMenu = XmCreatePulldownMenu(
	  menuBar,
	  "fileMenu",
	  NULL, 0);
  /* note that the fileMenu is NOT managed yet */

  /* create Button in menu that exits the application */
  quit = XtVaCreateManagedWidget(
	  "quit",
	  xmPushButtonWidgetClass,
	  fileMenu,
	  NULL);
  /* connect the menu with the filebutton */
  XtVaSetValues(fileButton,XmNsubMenuId,fileMenu,NULL);

  /* arrange for buttons to make their callbacks... */
  XtAddCallback(quit,XmNactivateCallback,Quit,0);
  XtAddCallback(resetPlotButton,XmNactivateCallback,ResetPlot,0);
  XtAddCallback(pauseButton,XmNactivateCallback,Pause,0);
  XtAddCallback(stopButton,XmNactivateCallback,Stop,0);
  XtAddCallback(showEquilButton,XmNactivateCallback,ShowEquil,0);

  XtAddCallback(platCommand,XmNcommandChangedCallback,CommandUpdate,0);
  XtAddCallback(platCommand,XmNcommandEnteredCallback,PlatCommand,0);

  /* Set Mainwindow areas */
  XmMainWindowSetAreas(mainWindow,menuBar,platCommand,NULL,NULL,workSpace);

  /* Set it all in motion... */
  XtAppAddActions(app_context, actions, XtNumber(actions));
  XtRealizeWidget(topLevel);

  XtAppMainLoop(app_context);
}

void
SetUpPixmap(Widget w)
{
  XGCValues values;

  bitmap.froth_bitmap = XCreatePixmap(
	XtDisplay(w),
	RootWindowOfScreen(XtScreen(w)),
	bitmap.pixmap_width,
	bitmap.pixmap_height, 1);
  values.foreground = 0;
  values.background = 1;
  bitmap.draw_gc = XCreateGC(
	XtDisplay(w),
	bitmap.froth_bitmap,
	GCForeground | GCBackground, &values);
  values.foreground = 1;
  values.background = 0;
  bitmap.undraw_gc = XCreateGC(
	XtDisplay(w),
	bitmap.froth_bitmap,
	GCForeground | GCBackground, &values);
  /* Pixmap must now be cleared... */
  XFillRectangle(
	XtDisplay(w),
	bitmap.froth_bitmap,
	bitmap.undraw_gc,
	0, 0,
	bitmap.pixmap_width+1,
	bitmap.pixmap_height+1);
}

void Test()
{
  Widget w;

  w = bitmap.frothPicture;
  XDrawLine(XtDisplay(w),
	bitmap.froth_bitmap, bitmap.draw_gc,
	100, 100, 200, 200);
  XDrawPoint(XtDisplay(w),
	bitmap.froth_bitmap, bitmap.draw_gc,
	20, 20);
/*
  XCopyPlane(XtDisplay(w),
	bitmap.froth_bitmap,
	XtWindow(w),
	DefaultGCOfScreen(XtScreen(w)), 0, 0,
	bitmap.pixmap_width,
	bitmap.pixmap_height,
	0, 0, 1);
*/
}


/*
 *
 * COMMAND DRIVER ROUTINES (PREVIOUSLY THE MAIN 'PLAT' ROUTINE):
 *
 */

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    P L A T _ C O M M A N D _ D R I V E R ( )    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	plat_command_driver()
 *
 *	Arguments:	none
 *
 *	Return value:	none
 *
 *	Action:		A plain (non-X window) front end to the program.
 *			This routine drives a simple command interface to
 *			the program (after doing a little bit of
 *			initialisation).  It is the interface which the user
 *			first encounters when he starts up the program.  The
 *			user can subsequently start up the X-window interface
 *			by giving the `g' command.  This `plain' interface is
 *			useful for running jobs in the background.
 *
 *			The details of the command interface can be learned
 *			by consulting `help.txt'.  Basically the commands are
 *			single letter commands.  Lower-case letters are simple
 *			commands and take no arguments.  Upper-case letters
 *			do take arguments.  A simple looping mechanism is also
 *			implemented e.g.
 *				{ <commands> }10
 *			would execute the commands in curly braces ten times.
 *			Loops may be nested.  This routine calls `sprompt()'
 *			to prompt for input and `s_execute()' to execute the
 *			commands in the input string.
 *
 *****************************************************************************/
void
plat_command_driver()
{
  void plerror(), setconstants(), trans(), vorspray();
  void mglinit(), mglline(), mglclose(), mglreset();
  void vfoamplot(), foamplot(), polyplot(), voronoi(), hexnet(),
       vortrans(), ptbox(), hencky(), s_execute();
  short i, j, k, nestlevel, slen;
  REAL x1, y1, x2, y2, d, eps;
  float eps0, f;
  char s[256], c;
  boolean sgood;
  boolean sprompt();

  do {
    nestlevel=0; slen=0; sgood=TRUE;
    do {
      if (sgood=sprompt(s+slen, nestlevel)) {
        slen=strlen(s);
        /* First check if the execution of the string ought to be suspended */
        /* pending further input */
        nestlevel=0;
        for (i=0; i<slen; i++) {
          switch(s[i]) {
          case '{' : nestlevel++; break;
          case '}' : nestlevel--;
                     if (nestlevel<0) {
                       plerror("too many right brackets");
                       sgood=FALSE; i=slen;
                     }
                     break;
          case '\\' :
            s[i]=' ';
            if (i!=slen) {
              plerror("illegal embedded backslash character \'\\\'");
              sgood=FALSE; i=slen;
            }
            break;
          }
        }
      }
      s[slen]='\n'; slen++; s[slen]='\0';
    } while (nestlevel>0 && sgood);
    if (sgood) s_execute(s,0);
    fflush(stdout);
  } while (! x_window_mode);
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *    S P R O M P T ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	boolean sprompt(char s[], short nestlevel)
 *
 *	Arguments:	s	= input string typed in by user
 *			nestlevel = the current level of nesting of curly
 *				braces `{ }' (which is the syntax of a simple
 *				`for' loop in the command interpreter).
 *
 *	Return value:	Returns false if error on reading input.
 *
 *	Action:		Prompts the user for input and returns a line of input
 *			in the string `s'.  The prompt also keeps track of the
 *			level of nesting which occurs due to pairs of braces
 *			`{ }' (which implement the `for' loop syntax).  If the
 *			user types an open brace which is not matched by a
 *			close brace on the same line then the `nestlevel' is
 *			incremented.  This means that execution of the
 *			commands is suspended until a matching close brace is
 *			encountered.  So long as the `nestlevel' is non-zero,
 *			the prompt will incorporate `nestlevel' open braces
 *			to indicate that it is waiting for the corresponding
 *			number of close braces.  E.g. a single open brace
 *			pending would give rise to a prompt `plat{>'.
 *
 *****************************************************************************/
boolean sprompt(s, nestlevel)
char *s;
short nestlevel;
{
  short i, len;
  void plerror();
  printf("plat");
  for (i=0; i<nestlevel; i++) printf("{");
  printf("> ");
  do {
    if (!gets(s)) {
      plerror("error while reading command input"); return FALSE;
    }
    for (len=strlen(s); len>0 && isspace(s[len-1]); len--);
    s[len]= '\0';
  } while (len!=0 && s[len-1]=='\\');
  return TRUE;
}

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *S _ E X E C U T E ( ) *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	s_execute(char s[], short nestlevel)
 *
 *	Arguments:	s	= input string to be executed
 *			nestlevel = the current level of nesting of curly
 *				braces `{ }' (which is the syntax of a simple
 *				`for' loop in the command interpreter).
 *
 *	Return value:	none
 *
 *	Action:		Parses and executes a command string `s'.  This
 *			subroutine also calls itself recursively in order to
 *			implement nested loops.  The commands have a simple
 *			syntax and consist of a single letter.  They are
 *			classed into Lower-case commands, which take no
 *			arguments, and Upper-case commands, which may take
 *			arguments.  White space between commands is optional
 *			unless it is necessary to remove ambiguity.  A full
 *			description of the command interface is given in
 *			the file `help.txt'.
 *
 *			The only slight refinement of the command interface
 *			is an ability to do simple `for' loops.  By putting
 *			commands inside curly braces e.g.
 *				{  <commands>  }<integer>
 *			the commands inside the curly braces will be executed
 *			<integer> times.  The loops may also be nested.
 *
 *****************************************************************************/
void s_execute(s, nestlevel)
char *s;
short nestlevel;
{
  void plerror(), setconstants(), trans(), vorspray();
  void mglinit(), mglline(), mglclose(), mglreset();
  void vfoamplot(), foamplot(), polyplot(), voronoi(), hexnet(),
       vortrans(), ptbox(), hencky(), s_execute(), diffuse(),
       makefraction(), setsconst(), foamin(), foamout(),
       infoopen(), infowrite(), infoclose(), setinfo(), constantout(),
       cadjm();
#ifdef IRIS
  void zoom(), mgltest();
#endif
  REAL equil(), phifn();
  boolean sprompt();
  short i, ii, j, k, len;
  int n1, n2, numlen;
  char c, tok[MSTRING];
  static char fname[50]="foam.pb";
  boolean sgood, flag, fflag;
  static REAL eps=0.01;
  static boolean infoflag=FALSE;
  float f;
  len=strlen(s);
  for (i=0; i<len; i++) {
    if (!isspace(c=s[i])) {
      if (c=='{') {
        ii=i+1;
        while (++i<len && s[i]!='}');
        if (i==len) {
          plerror("missing right bracket \'}\'"); return;
        }
        sscanf(s+i+1,"%d%n",&n1,&numlen);
	i += numlen;
        for (j=0; j<n1 && !xuser_stop; j++) s_execute(s+ii, nestlevel+1);
      }
      else if (c=='}') return;
      else {
        /* Execute simple commands */
        switch (c) {
        case 'p' : foamplot(0); break;
        case 'P' :
	  n1 = 1;
	  flag= FALSE;
	  mgl_hpgl_flag = mgl_tek_flag = mgl_ps_flag = FALSE;
          while (isspace(s[i+1]) && i<len) i++;
	  while (s[i+1]=='-' && i<len) {
	    i++;
	    if (strncmp(s+i+1,"n",1)==0) {
	      /* Do a multiple box plot... */
	      i++;
              while (isspace(s[i+1]) && i<len) i++;
	      if (!isdigit(s[i+1]) || s[i+1]=='-' ) {
		plerror("error: option -n must be followed by a positive integer");
		return;
	      }
              sscanf(s+i+1,"%d%n",&n1,&numlen);
	      i += numlen;
	    }
	    else if (strncmp(s+i+1,"hpgl",4)==0) {
	      i += 4;
              while (isspace(s[i+1]) && i<len) i++;
	      for (j=0; !isspace(s[i+1]) && i<len; i++, j++) {
		mgl_hpgl_filename[j] = s[i+1];
	      }
	      mgl_hpgl_filename[j] = '\0';
	      mgl_hpgl_flag = TRUE;
	      printf("writing file %s ...\n",mgl_hpgl_filename);
	    }
	    else if (strncmp(s+i+1,"tek",3)==0) {
	      i += 3;
              while (isspace(s[i+1]) && i<len) i++;
	      for (j=0; !isspace(s[i+1]) && i<len; i++, j++) {
		mgl_tek_filename[j] = s[i+1];
	      }
	      mgl_tek_filename[j] = '\0';
	      mgl_tek_flag = TRUE;
	      printf("writing file %s ...\n",mgl_tek_filename);
	    }
	    else if (strncmp(s+i+1,"ps",2)==0) {
	      i += 2;
              while (isspace(s[i+1]) && i<len) i++;
	      for (j=0; !isspace(s[i+1]) && i<len; i++, j++) {
		mgl_ps_filename[j] = s[i+1];
	      }
	      mgl_ps_filename[j] = '\0';
	      mgl_ps_flag = TRUE;
	      printf("writing file %s ...\n",mgl_ps_filename);
	    }
	    else if (strncmp(s+i+1,"hatch",5)==0) {
	      i += 5;
	      flag= TRUE;
	    }
	    else {
	      plerror("error: unrecognised plot option");
	      return;
	    }
            while (isspace(s[i+1]) && i<len) i++;
	  }
          polyplot((short) n1,flag); break;
	  /* Reset flags again */
	  mgl_hpgl_flag = mgl_tek_flag = mgl_ps_flag = FALSE;
        case 'a' : equil(1); break;
        case 'e' : j=0;
                   do {
         /*            if (xuser_show_equil==1) foamplot(0);
		     if (XtAppPending(app_context)) {
		       XtAppProcessEvent(app_context,XtIMXEvent);
		     }
		     while (xuser_pause && !xuser_stop) {
		       if (XtAppPending(app_context)) {
		         XtAppProcessEvent(app_context,XtIMXEvent);
		       }
		     };*/
		     j++;
                   } while (equil(0)>equilsup && ((REAL) j)<smaxiter
			    && !xuser_stop);
                   printf("network energy = %f\n",netenergy);
                   printf("%d iterations taken...\n",j);
                   break;
        case 'E' : j=0;
                   do {
                     if (xuser_show_equil != 0) foamplot(0);
		     if (XtAppPending(app_context)) {
		       XtAppProcessEvent(app_context,XtIMXEvent);
		     }
#ifdef IRIS
                     zoom();
#endif
		     while (xuser_pause && !xuser_stop) {
		       if (XtAppPending(app_context)) {
		         XtAppProcessEvent(app_context,XtIMXEvent);
		       }
		     };
		     j++;
                   } while (equil(0)>equilsup && ((REAL) j)<smaxiter
			    && !xuser_stop);
                   printf("network energy = %f\n",netenergy);
                   printf("%d iterations taken...\n",j);
                   break;
        case 'd' : diffuse(); break;
        case 'D' : diffuse(); break;
        case 'F' :
          while (isspace(s[i+1]) && i<len) i++;
          flag = fflag = FALSE;
          if (s[i+1]=='-' && isalpha(s[i+2])) {
            switch (s[i+2]) {
            case 'd' : fflag=TRUE; break;
            case 'e' : flag=FALSE; break;
            case 'E' : flag=TRUE; break;
            default : flag=FALSE; break;
            }
            i += 2;
          }
          sscanf(s+i+1,"%f%n",&f,&numlen);
	  i += numlen;
          if (fflag) {
            if (volfrac==0.0) {
              plerror("error: volume fraction not set");
              plerror("error: first use 'F <number> without the '-d' option");
              exit(0);
            }
            volfrac += (REAL)f;
          }
          else { volfrac= (REAL)f; }
          makefraction(volfrac,flag); break;
        case 'b' : hencky(eps); break;
        case 'B' :
          sscanf(s+i+1,"%f%n",&f,&numlen);
	  i += numlen;
          eps= (REAL) f;
          hencky(eps); break;
        case 'h' : case 'H' : system("more help.txt"); break;
        case 'g' : case 'G' :
#ifdef X_WINDOWS
	  if (! x_window_mode) x_window_mode = TRUE;
	  if (i < len-2) {
	    len = 0;
	    plerror("warning: rest of command string has been ignored");
	  }
#else
          ginteract=TRUE; mglreset();
#endif
          break;
        case 'i' :
          if (!infoflag) { infoopen(); infoflag=TRUE; }
          infowrite();
          break;
        case 'I' :
          flag=FALSE; fflag=TRUE;
          while (isspace(s[i+1]) && i<len) i++;
          if (s[i+1]=='-' && isalpha(s[i+2])) {
            switch (s[i+2]) {
            case 'w' : flag=TRUE; break;
            case 'W' : fflag=FALSE; break;
            case 'c' :
              if (infoflag) { infoclose(); infoflag=FALSE; }
              break;
            default : flag=FALSE; break;
            }
            i += 2;
          }
          if (s[i]!='c') {
            sscanf(s+i+1,"%s%n",tok,&numlen);
	    i += numlen;
            if (!infoflag) {
              setinfo(tok,flag,fflag);
            }
            else {
              plerror("cannot set new info while 'info.pb' is open");
              plerror("use 'I -c' option to close this file");
            }
          }
          break;
        case 'r' : foamin("foam.pb"); setconstants(); break;
        case 'R' :
          sscanf(s+i+1,"%s%n",fname,&numlen);
	  i += numlen;
          foamin(fname); setconstants();
          break;
        case 'w' : foamout("foam.pb"); break;
        case 'W' :
          sscanf(s+i+1,"%s%n",fname,&numlen);
	  i += numlen;
          foamout(fname);
          break;
        case 'u' : case 'U' : break;
        case 'n' : mglreset(); hexnet(1,1,0.4); setconstants(); break;
        case 'N' :
          mglreset();
          sscanf(s+i+1,"%d%d%f%n",&n1,&n2,&f,&numlen);
	  i += numlen;
          hexnet((short) n1,(short) n2,(REAL) f); setconstants(); break;
        case 'v' : 
          mglreset(); nc=MCELL; voronoi(nc); setconstants(); break;
        case 'V' :
          mglreset();
          while (isspace(s[i+1]) && i<len) i++;
          if (s[i+1]=='-' && isalpha(s[i+2])) {
            switch (s[i+2]) {
            case 'r' :
              i += 2;
              sscanf(s+i+1,"%f%n",&f,&numlen);
	      i += numlen;
	      srfrac=(REAL) f;
              break;
            default : printf("unrecognised option\n"); i += 2; break;
            }
          }
          sscanf(s+i+1,"%d%n",&n1,&numlen);
	  i += numlen;
          nc=(short) n1; voronoi(nc); setconstants(); break;
        case 's' : break;  /* not yet implemented */
        case 'S' :
          while (isspace(s[i+1]) && i<len) i++;
          if (s[i+1]=='-' && isalpha(s[i+2])) {
            switch (s[i+2]) {
            case 'w' :
              i += 2;
              sscanf(s+i+1,"%s%n",tok,&numlen);
              i += numlen;
              constantout(tok);
              break;
            default : printf("unrecognised option\n"); i += 2; break;
            }
          }
          else {
            sscanf(s+i+1,"%s%f%n",tok,&f,&numlen);
            i += numlen;
            setsconst(tok,(REAL) f);
          }
          break;
        case 't' : break; /* Not implemented */
        case 'q' : case 'Q' :
          if (infoflag) infoclose();
          plerror("gurgle, gurgle, gurgle...");
          exit(0); break;
        case 'x' : case 'X' :
          if (nestlevel==0) {
            if (infoflag) infoclose();
            plerror("gurgle, gurgle, gurgle..."); exit(0);
          }
          else return;
        case 'z' :
#if defined(IRIS)
          zoom(); foamplot(0);
#else
          printf("Iris: zoom not implemented\n");
#endif
          break;
        case 'Z' :
#if defined(IRIS)
          zoom(); foamplot(0);
#else
          printf("Iris: zoom not implemented\n");
#endif
          break;
        default : ;
          /*plerror("unrecognised command letter");*/
        }
      }
    }
  }
}
