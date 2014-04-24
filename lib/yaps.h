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
#include <stdarg.h>

void yaps_message(const char *fmt, ...);
void yaps_quit(const char *fmt, ...);
void yaps_sysquit(const char *fmt, ...);
void yaps_yapper(void (*yapper)(const char *format, va_list ap));
