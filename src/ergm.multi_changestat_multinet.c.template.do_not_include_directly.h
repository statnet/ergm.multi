/*  File src/wtchangestats_multinet.c in package ergm.multi, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */

#include "ergm_storage.h"


ETYPE(I_CHANGESTAT_FN)(ETYPE(i__, subnets)){
  int *iinputs = IINPUT_PARAM;
  ALLOC_AUX_STORAGE(1, ETYPE(Store, Subnets), sn);
  sn->ns = *(iinputs++);
  sn->inwp = nwp;
  sn->onwp = R_Calloc(sn->ns, ETYPE(Network) *);
  sn->onwp--; // The -- is because ETYPE(Network) IDs count from 1.

  /* Set up the layer information. */
  sn->sid = (Vertex *) iinputs - 1; // The -1 is because Vertex IDs count from 1.
  iinputs += N_NODES;
  sn->smap = (Vertex *) iinputs - 1;
  iinputs += N_NODES;

  for(unsigned int i=1; i<=sn->ns; i++){
    Vertex lnnodes, lbip;
    if(BIPARTITE){
      lbip = lnnodes = *(iinputs++);
      lnnodes += *(iinputs++);
    }else{
      lbip = 0;
      lnnodes = *(iinputs++);
    }

    sn->onwp[i] = ETYPE(NetworkInitialize)(NULL, NULL, IFEWT(NULL,) 0, lnnodes, DIRECTED, lbip);
  }

  IFELSEEWT(
            EXEC_THROUGH_NET_EDGES_PRE(t, h, e, w, {
               ETYPE(SetEdge)(MN_IO_TAIL(sn, t), MN_IO_HEAD(sn, h), w, sn->onwp[MN_SID_TAIL(sn, t)]);
             }),
           EXEC_THROUGH_NET_EDGES_PRE(t, e, h, {
               ToggleKnownEdge(MN_IO_TAIL(sn, t), MN_IO_HEAD(sn, h),sn->onwp[MN_SID_TAIL(sn, t)], FALSE);
             })
           );
}

ETYPE(U_CHANGESTAT_FN)(ETYPE(u__, subnets)){ 
  GET_AUX_STORAGE(ETYPE(Store, Subnets), sn);
  IFELSEEWT(
            ETYPE(SetEdge)(MN_IO_TAIL(sn, tail), MN_IO_HEAD(sn, head), weight, sn->onwp[MN_SID_TAIL(sn, tail)]),
            ToggleKnownEdge(MN_IO_TAIL(sn, tail), MN_IO_HEAD(sn, head),sn->onwp[MN_SID_TAIL(sn, tail)], edgestate)
            );
}

ETYPE(F_CHANGESTAT_FN)(ETYPE(f__, subnets)){
  GET_AUX_STORAGE(ETYPE(Store, Subnets), sn);
  for(unsigned int i=1; i<=sn->ns; i++)
    ETYPE(NetworkDestroy)(sn->onwp[i]);
  sn->onwp++;
  R_Free(sn->onwp);
}

// MultiNet: Take a weighted networkwise sum of the networks' statistics.

ETYPE(I_CHANGESTAT_FN)(ETYPE(i_, MultiNet)){
  /*
    iinputs expects:
    1: number of weights (nwts)
    inputs expects:
    nwts*ns: matrix of weights, in network-major order
  */
  
  GET_AUX_STORAGE(ETYPE(Store, Subnets), sn);
  unsigned int ns = sn->ns;
  unsigned int nwts = *IINPUT_PARAM;
  double *wts = INPUT_PARAM;
  
  ALLOC_STORAGE(ns, ETYPE(Model)*, ms);

  SEXP submodels = getListElement(mtp->R, "submodels");
  unsigned int submodpos = 0;
  for(unsigned int i=1; i<=sn->ns; i++){
    unsigned int used=FALSE;
    for(unsigned int j=0; j<nwts; j++){
      if(wts[j]!=0){
	used=TRUE;
	break;
      }
    }
    wts += nwts; // OK to clobber it here.
    if(used){
      ms[i-1] = ETYPE(ModelInitialize)(VECTOR_ELT(submodels, submodpos), NULL, sn->onwp[i], FALSE);
      submodpos++;
    }else ms[i-1] = NULL;
  }
  ETYPE(DELETE_IF_UNUSED_IN_SUBMODELS)(u_func, ms, sn->ns);
  ETYPE(DELETE_IF_UNUSED_IN_SUBMODELS)(z_func, ms, sn->ns);
}

ETYPE(C_CHANGESTAT_FN)(ETYPE(c_, MultiNet)){
  GET_AUX_STORAGE(ETYPE(Store, Subnets), sn);
  GET_STORAGE(ETYPE(Model)*, ms);
  unsigned int nwts = *IINPUT_PARAM;
  double *wts = INPUT_PARAM;

  unsigned int i = MN_SID_TAIL(sn, tail);
  ETYPE(Model) *m = ms[i-1];
  if(m){ // NULL if network has weights 0.
    Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
    ETYPE(ChangeStats1)(st, sh, IFEWT(weight,) sn->onwp[i], m, edgestate);

    wts += (i-1)*nwts; // Position of that network's weight vector.
    for(unsigned int j=0; j<m->n_stats; j++)
      for(unsigned int k=0; k<nwts; k++)
	CHANGE_STAT[j*nwts+k] += m->workspace[j]*wts[k];
  }
}

ETYPE(Z_CHANGESTAT_FN)(ETYPE(z_, MultiNet)){
  GET_AUX_STORAGE(ETYPE(Store, Subnets), sn);
  GET_STORAGE(ETYPE(Model)*, ms);
  unsigned int nwts = *IINPUT_PARAM;
  double *wts = INPUT_PARAM;

  for(unsigned int i=1; i<=sn->ns; i++){
    ETYPE(Model) *m = ms[i-1];
    if(m){ // NULL if network has weights 0.
      ETYPE(ZStats)(sn->onwp[i], m, FALSE);

      wts += (i-1)*nwts; // Position of that network's weight vector.
      for(unsigned int j=0; j<m->n_stats; j++)
        for(unsigned int k=0; k<nwts; k++)
          CHANGE_STAT[j*nwts+k] += m->workspace[j]*wts[k];
    }
  }
}

ETYPE(F_CHANGESTAT_FN)(ETYPE(f_, MultiNet)){
  GET_AUX_STORAGE(ETYPE(Store, Subnets), sn);
  GET_STORAGE(ETYPE(Model)*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    if(ms[i-1]) ETYPE(ModelDestroy)(sn->onwp[i], ms[i-1]);
  }
}

// MultiNets: Concatenate the networks' statistics; network statistic counts may be heterogeneous.

ETYPE(I_CHANGESTAT_FN)(ETYPE(i_, MultiNets)){
  int *iinputs = IINPUT_PARAM; 
  GET_AUX_STORAGE(ETYPE(Store, Subnets), sn);
  unsigned int ns = sn->ns;
  unsigned int *pos = (unsigned int *) iinputs;
  ALLOC_STORAGE(ns, ETYPE(Model)*, ms);

  SEXP submodels = getListElement(mtp->R, "submodels");
  unsigned int submodpos = 0;
  for(unsigned int i=1; i<=sn->ns; i++){
    if(pos[i-1]!=pos[i]){
      ms[i-1] = ETYPE(ModelInitialize)(VECTOR_ELT(submodels, submodpos), NULL, sn->onwp[i], FALSE);
      submodpos++;
    }
  }
  ETYPE(DELETE_IF_UNUSED_IN_SUBMODELS)(u_func, ms, sn->ns);
  ETYPE(DELETE_IF_UNUSED_IN_SUBMODELS)(z_func, ms, sn->ns);
}

ETYPE(C_CHANGESTAT_FN)(ETYPE(c_, MultiNets)){
  unsigned int *pos = (unsigned int *) IINPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(ETYPE(Store, Subnets), sn);
  GET_STORAGE(ETYPE(Model)*, ms);

  unsigned int i = MN_SID_TAIL(sn, tail);
  Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
  if(pos[i-1]!=pos[i]){
    ETYPE(Model) *m = ms[i-1];
    ETYPE(ChangeStats1)(st, sh, IFEWT(weight,) sn->onwp[i], m, edgestate);
    memcpy(CHANGE_STAT + (unsigned int)(pos[i-1]), m->workspace, m->n_stats*sizeof(double));
  }
}

ETYPE(Z_CHANGESTAT_FN)(ETYPE(z_, MultiNets)){
  unsigned int *pos = (unsigned int *) IINPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(ETYPE(Store, Subnets), sn);
  GET_STORAGE(ETYPE(Model)*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    if(pos[i-1]!=pos[i]){
      ETYPE(Model) *m = ms[i-1];
      ETYPE(ZStats)(sn->onwp[i], m, FALSE);
      memcpy(CHANGE_STAT + (unsigned int)(pos[i-1]), m->workspace, m->n_stats*sizeof(double));
    }
  }
}

ETYPE(F_CHANGESTAT_FN)(ETYPE(f_, MultiNets)){
  unsigned int *pos = (unsigned int *) IINPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(ETYPE(Store, Subnets), sn);
  GET_STORAGE(ETYPE(Model)*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    if(pos[i-1]!=pos[i]){
      ETYPE(ModelDestroy)(sn->onwp[i], ms[i-1]);
    }
  }
}

// wtByNetDStats

ETYPE(I_CHANGESTAT_FN)(ETYPE(i_, ByNetDStats)){
  ETYPE(Model) *m = STORAGE = ETYPE(ModelInitialize)(getListElement(mtp->R, "submodel"), NULL, nwp, FALSE);
  ETYPE(DELETE_IF_UNUSED_IN_SUBMODEL)(u_func, m);
  ETYPE(DELETE_IF_UNUSED_IN_SUBMODEL)(z_func, m);
}

ETYPE(C_CHANGESTAT_FN)(ETYPE(c_, ByNetDStats)){
  unsigned int *pos = (unsigned int *) IINPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(ETYPE(Store, Subnets), sn);
  GET_STORAGE(ETYPE(Model), m);

  unsigned int i = MN_SID_TAIL(sn, tail);
  if(pos[i-1]!=pos[i]){
    ETYPE(ChangeStats1)(tail, head, IFEWT(weight,) nwp, m, edgestate);
    memcpy(CHANGE_STAT + (unsigned int)pos[i], m->workspace, m->n_stats*sizeof(double));
  }
}

ETYPE(Z_CHANGESTAT_FN)(ETYPE(z_, ByNetDStats)){
  unsigned int *pos = (unsigned int *) IINPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(ETYPE(Store, Subnets), sn);
  GET_STORAGE(ETYPE(Model), m);

  for(unsigned int i=1; i<=sn->ns; i++)
    if(pos[i-1]!=pos[i]){
      ETYPE(ZStats)(nwp, m, FALSE);
      memcpy(CHANGE_STAT + (unsigned int)pos[i], m->workspace, m->n_stats*sizeof(double));
    }
}

ETYPE(F_CHANGESTAT_FN)(ETYPE(f_, ByNetDStats)){
  GET_STORAGE(ETYPE(Model), m);
  ETYPE(ModelDestroy)(nwp, m);
  STORAGE = NULL;
}
