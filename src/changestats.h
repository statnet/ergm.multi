#ifndef _CHANGESTATS_H_
#define _CHANGESTATS_H_

#include "ergm_Rutil.h"
#include "ergm_changestat.h"

static inline void cutoff_error(ModelTerm *mtp){
  error("%s", CHAR(STRING_ELT(getListElement(mtp->R, "cutoff.message"), 0)));
}

#endif // _CHANGESTATS_H_
