/*  File src/ergm.multi_changestat_multinet.h in package ergm.multi, part of
 *  the Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _ERGM_MULTI_CHANGESTAT_MULTINET_H_
#define _ERGM_MULTI_CHANGESTAT_MULTINET_H_

#include "ergm_edgetree.h"

#include "ergm_edgetype_set_binary.h"
#include "ergm.multi_changestat_multinet.h.template.do_not_include_directly.h"

#define MN_IO_TAIL(sn, t) ((Vertex) ((sn)->smap[t]))
#define MN_IO_HEAD(sn, h) ((Vertex) ((sn)->smap[h]))
#define MN_SID_TAIL(sn, t) ((Vertex) ((sn)->sid[t]))
#define MN_SID_HEAD(sn, h) ((Vertex) ((sn)->sid[h]))

#define MN_IGETWT(sn, l,a,b) (GetEdge(MN_OI_TAIL(sn, l, a), MN_OI_HEAD(sn, l, b), sn->inwp))
#define MN_ISETWT(sn, l,a,b,w) (SetEdge(MN_OI_TAIL(sn, l, a), MN_OI_HEAD(sn, l, b),w,(sn)->inwp))

#endif // _ERGM_MULTI_CHANGESTAT_MULTINET_H_
