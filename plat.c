#include "include.h"
#include <string.h>
#include <ctype.h>

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
      dvy[MVERT], dcp[MCELL];
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
     "nbar", "edgefrac", "z", "arootbar", "minenergy", "adjmat"};
boolean info_list[MINFO], info_wlist[MINFO], foamlike, ginteract=FALSE,
        notopol=FALSE;

/*****************************************************************************
 ***     * *    * *      * *                          * *    * *    * *    ***
 ***    *   * *    *  *     *        M A I N ( )   *     * *    * *    *   ***
 ***   *    *      *        *                   *        *      *      *   ***
 *****************************************************************************
 *
 *	Subroutine:	main()
 *
 *	Action:		The main routine drives a simple command interface to
 *			the program (after doing a little bit of
 *			initialisation).  If I were writing it now I would
 *			use `Tcl' as it is designed for exactly this sort of
 *			thing!
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
main()
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
  mglreset();
  nc=MCELL;
  for (i=0; i<MINFO; i++) {
    info_list[i]=FALSE; info_wlist[i]=FALSE;
  }
  for (;;) {
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
  }
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
       infoopen(), infowrite(), infoclose(), setinfo(), constantout();
#ifdef IRIS
  void zoom(), mgltest();
#endif
  REAL equil(), phifn();
  boolean sprompt();
  short i, ii, j, k, len;
  int n1, n2;
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
        sscanf(s+i+1,"%d",&n1);
        while (!isalpha(s[i+1]) && i<len) i++;
        for (j=0; j<n1; j++) s_execute(s+ii, nestlevel+1);
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
              sscanf(s+i+1,"%d",&n1);
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
            while (!isalpha(s[i+1]) && s[i+1]!='-' && i<len) i++;
	  }
	  /* oldstuff */
          polyplot((short) n1,flag); break;
	  /* Reset flags again */
	  mgl_hpgl_flag = mgl_tek_flag = mgl_ps_flag = FALSE;
        case 'a' : equil(1); break;
        case 'e' : j=0;
                   do { j++;
                   } while (equil(0)>equilsup && ((REAL) j)<smaxiter);
                   printf("network energy = %f\n",netenergy);
                   printf("%d iterations taken...\n",j);
                   break;
        case 'E' : j=0;
                   do {
#ifdef IRIS
                     zoom();
#endif
                     foamplot(0); j++;
                   } while (equil(0)>equilsup && ((REAL) j)<smaxiter);
                   printf("network energy = %f\n",netenergy);
                   printf("%d iterations taken...\n",j);
                   break;
        case 'd' : diffuse(); break;
        case 'D' : diffuse(); break;
        case 'F' :
          while (isspace(s[i+1]) && i<len) i++;
          fflag=FALSE;
          if (s[i+1]=='-' && isalpha(s[i+2])) {
            switch (s[i+2]) {
            case 'd' : fflag=TRUE; break;
            case 'e' : flag=FALSE; break;
            case 'E' : flag=TRUE; break;
            default : flag=FALSE; break;
            }
            i += 2;
          }
          sscanf(s+i+1,"%f",&f);
          if (fflag) {
            if (volfrac==0.0) {
              plerror("error: volume fraction not set");
              plerror("error: first use 'F <number> without the '-d' option");
              exit(0);
            }
            volfrac += (REAL)f;
          }
          else { volfrac= (REAL)f; }
          while (!isalpha(s[i+1]) && i<len) i++;
          makefraction(volfrac,flag); break;
        case 'b' : hencky(eps); break;
        case 'B' :
          sscanf(s+i+1,"%f",&f);
          while (!isalpha(s[i+1]) && i<len) i++;
          eps= (REAL) f;
          hencky(eps); break;
        case 'h' : case 'H' : system("more help.txt"); break;
        case 'g' : case 'G' :
          ginteract=TRUE; mglreset();
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
            sscanf(s+i+1,"%s",tok);
            while (isspace(s[i+1]) && i<len) i++;
            i += strlen(tok);
            while (!isalpha(s[i+1]) && i<len) i++;
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
          sscanf(s+i+1,"%s",fname);
          while (isspace(s[i+1]) && i<len) i++;
          i += strlen(fname);
          while (!isalpha(s[i+1]) && i<len) i++;
          foamin(fname); setconstants();
          break;
        case 'w' : foamout("foam.pb"); break;
        case 'W' :
          sscanf(s+i+1,"%s",fname);
          while (isspace(s[i+1]) && i<len) i++;
          i += strlen(fname);
          while (!isalpha(s[i+1]) && i<len) i++;
          foamout(fname);
          break;
        case 'u' : case 'U' : break;
        case 'n' : mglreset(); hexnet(1,1,0.4); setconstants(); break;
        case 'N' :
          mglreset();
          sscanf(s+i+1,"%d%d%f",&n1,&n2,&f);
          while (!isalpha(s[i+1]) && i<len) i++;
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
              while (isspace(s[i+1]) && i<len) i++;
              sscanf(s+i+1,"%f",&f); srfrac=(REAL) f;
              while (!isspace(s[i+1]) && i<len) i++;
              while (isspace(s[i+1]) && i<len) i++;
              break;
            default : printf("unrecognised option\n"); i += 2; break;
            }
          }
          sscanf(s+i+1,"%d",&n1);
          while (!isalpha(s[i+1]) && i<len) i++;
          nc=(short) n1; voronoi(nc); setconstants(); break;
        case 's' : break;  /* not yet implemented */
        case 'S' :
          while (isspace(s[i+1]) && i<len) i++;
          if (s[i+1]=='-' && isalpha(s[i+2])) {
            switch (s[i+2]) {
            case 'w' :
              i += 2;
              while (isspace(s[i+1]) && i<len) i++;
              sscanf(s+i+1,"%s",tok);
              i += strlen(tok);
              constantout(tok);
              break;
            default : printf("unrecognised option\n"); i += 2; break;
            }
          }
          else {
            sscanf(s+i+1,"%s%f",tok,&f);
            while (isspace(s[i+1]) && i<len) i++;
            i += strlen(tok);
            while (!isalpha(s[i+1]) && i<len) i++;
            setsconst(tok,(REAL) f);
          }
          break;
        case 't' :
#ifdef IRIS
          mgltest();
#endif
          break;
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
          printf("zoom not implemented\n");
#endif
          break;
        case 'Z' :
#if defined(IRIS)
          zoom(); foamplot(0);
#else
          printf("zoom not implemented\n");
#endif
          break;
        default : ;
          /*plerror("unrecognised command letter");*/
        }
      }
    }
  }
}
