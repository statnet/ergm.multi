/*  File src/ergm.multi_wtchangestat_multinet.h in package ergm.multi, part of
 *  the Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _ERGM_MULTI_WTCHANGESTAT_MULTINET_H_
#define _ERGM_MULTI_WTCHANGESTAT_MULTINET_H_

#include "ergm_wtedgetree.h"

#include "ergm_edgetype_set_double.h"
#include "ergm.multi_changestat_multinet.h.template.do_not_include_directly.h"

#define WtMN_IO_TAIL(sn, t) ((Vertex) ((sn)->smap[t]))
#define WtMN_IO_HEAD(sn, h) ((Vertex) ((sn)->smap[h]))
#define WtMN_SID_TAIL(sn, t) ((Vertex) ((sn)->sid[t]))
#define WtMN_SID_HEAD(sn, h) ((Vertex) ((sn)->sid[h]))

#define WtMN_IGETWT(sn, l,a,b) (WtGetEdge(WtMN_OI_TAIL(sn, l, a), WtMN_OI_HEAD(sn, l, b), sn->inwp))
#define WtMN_ISETWT(sn, l,a,b,w) (WtSetEdge(WtMN_OI_TAIL(sn, l, a), WtMN_OI_HEAD(sn, l, b),w,(sn)->inwp))

/* If STRICT_Wt_HEADERS is not set, give the terms more generic names. */
#ifndef STRICT_Wt_HEADERS

#define MN_IO_TAIL WtMN_IO_TAIL
#define MN_IO_HEAD WtMN_IO_HEAD
#define MN_SID_TAIL WtMN_SID_TAIL
#define MN_SID_HEAD WtMN_SID_HEAD

#define MN_IGETWT WtMN_IGETWT
#define MN_ISETWT WtMN_ISETWT

#endif // STRICT_Wt_HEADERS

#endif // _ERGM_MULTI_WTCHANGESTAT_MULTINET_H_
