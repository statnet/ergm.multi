/*  File src/ergm.multi_wtchangestat_multinet.h in package ergm.multi, part of
 *  the Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "ergm.multi_changestat_multinet_common.do_not_include_directly.h"

/* network-aware macros eponymous to ergm_changestat.h */
/* These routines need a way to specify which of the subnetworks is
   being affected, so they are disabled for now. */
/* #define MN_IS_OUTEDGE(sn, a,b) (EdgetreeSearch((a),(b),(sn)->onwp->outedges)!=0?1:0) */
/* #define MN_IS_INEDGE(sn, a,b) (EdgetreeSearch((a),(b),(sn)->onwp->inedges)!=0?1:0) */
/* #define MN_IS_UNDIRECTED_EDGE(sn, a,b) MN_IS_OUTEDGE((sn), MIN(a,b), MAX(a,b)) */
/* #define MN_MIN_OUTEDGE(sn, a) (EdgetreeMinimum((sn)->onwp->outedges, (a))) */
/* #define MN_MIN_INEDGE(sn, a) (EdgetreeMinimum((sn)->onwp->inedges, (a))) */
/* #define MN_NEXT_OUTEDGE(sn, e) (EdgetreeSuccessor((sn)->onwp->outedges,(e))) */
/* #define MN_NEXT_INEDGE(sn, e) (EdgetreeSuccessor((sn)->onwp->inedges,(e))) */
/* #define MN_STEP_THROUGH_OUTEDGES(sn, a,e,v) for((e)=MN_MIN_OUTEDGE((sn), a);((v)=MN_OUTVAL((sn), e))!=0;(e)=MN_NEXT_OUTEDGE((sn), e)) */
/* #define MN_STEP_THROUGH_INEDGES(sn, a,e,v) for((e)=MN_MIN_INEDGE((sn), a);((v)=MN_INVAL((sn), e))!=0;(e)=MN_NEXT_INEDGE((sn), e)) */
/* #define MN_STEP_THROUGH_OUTEDGES_DECL(sn, a,e,v) for(Edge e=MN_MIN_OUTEDGE((sn), a);MN_OUTVAL((sn), e)!=0;e=MN_NEXT_OUTEDGE((sn), e)) */
/* #define MN_STEP_THROUGH_INEDGES_DECL(sn, a,e,v) for(Edge e=MN_MIN_INEDGE((sn), a);MN_INVAL((sn), e)!=0;e=MN_NEXT_INEDGE((sn), e)) */
/* #define MN_EXEC_THROUGH_OUTEDGES(sn, a,e,v,subroutine) if(MN_DIRECTED((sn))){ MN_STEP_THROUGH_OUTEDGES_DECL((sn), a,e,v) {Vertex v=MN_OUTVAL((sn), e); subroutine} } else { MN_EXEC_THROUGH_EDGES((sn), a,e,v,subroutine) } */
/* #define MN_EXEC_THROUGH_INEDGES(sn, a,e,v,subroutine) if(MN_DIRECTED((sn))){ MN_STEP_THROUGH_INEDGES_DECL((sn), a,e,v) {Vertex v=MN_INVAL((sn), e); subroutine} } else { MN_EXEC_THROUGH_EDGES((sn), a,e,v,subroutine) } */
/* #define MN_EXEC_THROUGH_EDGES(sn, a,e,v,subroutine) { MN_STEP_THROUGH_OUTEDGES_DECL((sn), a,e,v) {Vertex v=MN_OUTVAL((sn), e); subroutine}  MN_STEP_THROUGH_INEDGES_DECL((sn), a,e,v) {Vertex v=MN_INVAL((sn), e); subroutine} } */
/* #define MN_EXEC_THROUGH_FOUTEDGES(sn, a,e,v,subroutine) MN_STEP_THROUGH_OUTEDGES_DECL((sn), a,e,v) {Vertex v=MN_OUTVAL((sn), e); subroutine} */
/* #define MN_EXEC_THROUGH_FINEDGES(sn, a,e,v,subroutine) MN_STEP_THROUGH_INEDGES_DECL((sn), a,e,v) {Vertex v=MN_INVAL((sn), e); subroutine} */
/* #define MN_EXEC_THROUGH_NET_EDGES(sn, a,b,e,subroutine) for(Vertex a=1; a <= N_NODES; a++)  MN_EXEC_THROUGH_FOUTEDGES((sn), a, e, b, {subroutine}); */
/* #define MN_TOGGLE(sn, a,b) (ToggleEdge((a),(b),(sn)->onwp)); */
/* #define MN_TOGGLE_DISCORD(sn, a,b) (ToggleEdge((a),(b),(sn)->onwp+1)); */
/* #define MN_GETWT(sn, a,b) (GetEdge(a,b,(sn)->onwp)) */
/* #define MN_SETWT(sn, a,b,w) (SetEdge(a,b,w,(sn)->onwp)) */

typedef struct {
  unsigned int ns;
  ETYPE(Network) *inwp, **onwp;
  Vertex *sid;
  Vertex *smap;
} ETYPE(Store, Subnets);
