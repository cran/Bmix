/********************************************************************************
 *
 * Bayesian Regression and Adaptive Sampling with Gaussian Process Trees
 * Copyright (C) 2005, University of California
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Questions? Contact Robert B. Gramacy (rbgramacy@ams.ucsc.edu)
 *
 ********************************************************************************/

#include "rhelp.h"
#ifdef RPRINT
#include <R_ext/Print.h>
FILE *bobbys_stdout = (FILE*) 0;
FILE *bobbys_stderr = (FILE*) 1;
#else
FILE *bobbys_stdout = stdio;
FILE *bobbys_stderr = stderr;
#endif
#include <stdarg.h>
#include <time.h>
#include <assert.h>

/*
 * bobbys_printf:
 *
 * a function many different types of printing-- in particular, using
 * the Rprintf if the code happens to be compiled with RPRINT,
 * othersie fprintf (takes the same arguments as fprintf)
 */

void bobbys_printf(FILE *outfile, const char *str, ...)
{
  va_list argp;
  va_start(argp, str);

  #ifdef RPRINT
  if(outfile == bobbys_stdout) Rvprintf(str, argp);
  else if(outfile == bobbys_stderr) REvprintf(str, argp);
  else vfprintf(outfile, str, argp);
  #else
  vfprintf(outfile, str, argp);
  #endif

  va_end(argp);
}


#ifndef RPRINT
/*
 * error:
 *
 * printf style function that reports errors to stderr
 */

void error(const char *str, ...)
{
  va_list argp;
  va_start(argp, str);

  bobbys_printf(stderr, "ERROR: ");
  vfprintf(stderr, str, argp);

  va_end(argp);
  bobbys_flush(stderr);

  /* add a final newline */
  bobbys_printf(stderr, "\n");

  /* kill the code */
  assert(0);
}


/*
 * warning:
 *
 * printf style function that reports warnings to stderr
 */

void warning(const char *str, ...)
{
  va_list argp;
  va_start(argp, str);

  bobbys_printf(stderr, "WARNING: ");
  vfprintf(stderr, str, argp);

  va_end(argp);
  bobbys_flush(stderr);

  /* add a final newline */
  bobbys_printf(stderr, "\n");
}
#endif

/*
 * bobbys_flush:
 *
 * a function for many different types of flushing--  in particular,
 * using * the R_FlushConsole the code happens to be compiled with
 * RPRINT, otherwise fflush
 */

void bobbys_flush(FILE *outfile)
{
#ifdef RPRINT
  R_FlushConsole();
#else
  fflush(outfile);
#endif
}


/*
 * bobbys__r_process_events:
 *
 * at least every 1 second(s) pass control back to
 * R so that it can check for interrupts and/or
 * process other R-gui events
 */

time_t bobbys__r_process_events(time_t itime)
{
#ifdef RPRINT
  time_t ntime = time(NULL);

  if(ntime - itime > 1) {
    R_FlushConsole();
    R_CheckUserInterrupt();
#if  (defined(HAVE_AQUA) || defined(Win32) || defined(Win64))
    R_ProcessEvents();
#endif
    itime = ntime;
  }
#endif
  return itime;
}

/* assert that uses proper R errors under RPRINT */

void bobbys_assert(int cnd){
  if(!cnd)
	error("failed assertion."); }
