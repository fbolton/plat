/* The Original of 'machine.h' has been lost in transit.  This is then */
/* (hopefully) a correct replacement.                                  */

/* Select Double Precision... */
#define DOUBLEPRECISION

#ifdef DOUBLEPRECISION
#define REAL double
#else
#define REAL float
#endif

/* Include the following line to generate a file of Hewlett-Packard Graphics
Language.  This is written to the file 'hpgl.p' */
/* #define HPGL */

/* Include the following line to generate a file for use with the O'Reilly
Hewlett-Packard Printer */
/* #define ORIHPGL */

/* Include the following line to use an Iris Workstation */
/* #define IRIS */

/* Include the following line to generate a Postscript file */
#define POSTSCRIPT

/* Include the following line to use 'X' windows */
#define X_WINDOWS
