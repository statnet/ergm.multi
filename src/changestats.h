/*  File src/changestats.h in package ergm.multi, part of the Statnet suite of
 *  packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _CHANGESTATS_H_
#define _CHANGESTATS_H_

#include "ergm_Rutil.h"
#include "ergm_changestat.h"

static inline void cutoff_error(ModelTerm *mtp){
  error("%s", CHAR(STRING_ELT(getListElement(mtp->R, "cutoff.message"), 0)));
}

#endif // _CHANGESTATS_H_
