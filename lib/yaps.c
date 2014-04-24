/*
 * Error handling
 * Copyright (C) 2011 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *     
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "yaps.h"

/*
 *    setup so user can define their own yapper
 */
static void (*yaps)(const char *format, va_list ap) = NULL;
static void call_yaps(const char *fmt, ...) {
  va_list ap;
  if ( yaps==NULL )
    yaps_quit("Cannot call_yaps on NULL yapper\n");
  va_start(ap, fmt);
  yaps(fmt, ap);
  va_end(ap);
}
void yaps_yapper(void (*yapper)(const char *format, va_list ap)) {
  yaps = yapper;
}

/*
 * Nonfatal message 
 */
void yaps_message(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  if ( yaps ) 
    yaps(fmt, ap);
  else
    vfprintf(stderr, fmt, ap);
  va_end(ap);
}
/*
 * Fatal message 
 */
void yaps_quit(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  if ( yaps ) 
    yaps(fmt, ap);
  else
    vfprintf(stderr, fmt, ap);
  va_end(ap);
  exit(1);
}
/*
 * Fatal error related to a system call.
 */
void yaps_sysquit(const char *fmt, ...)
{
  va_list ap;
  if ( yaps ) 
    call_yaps("%s: ", strerror(errno));
  else
    fprintf(stderr, "%s: ", strerror(errno));
  va_start(ap, fmt);
  if ( yaps ) 
    yaps(fmt, ap);
  else
    vfprintf(stderr, fmt, ap);
  va_end(ap);
  exit(1);
}
