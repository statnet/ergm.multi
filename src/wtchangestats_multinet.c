#include "ergm.multi_wtchangestat_multinet.h"
#include "ergm_wtchangestat_operator.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtmodel.h"
#include "ergm_storage.h"


WtI_CHANGESTAT_FN(i__wtsubnets){
  int *iinputs = IINPUT_PARAM;
  ALLOC_AUX_STORAGE(1, StoreWtSubnets, sn);
  sn->ns = *(iinputs++);
  sn->inwp = nwp;
  sn->onwp = Calloc(sn->ns, WtNetwork *);
  sn->onwp--; // The -- is because WtNetwork IDs count from 1.

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

    sn->onwp[i] = WtNetworkInitialize(NULL, NULL, NULL, 0, lnnodes, DIRECTED, lbip, 0, 0, NULL);
  }
  
  EXEC_THROUGH_NET_EDGES_PRE(t, h, e, w, {
      WtSetEdge(MN_IO_TAIL(sn, t), MN_IO_HEAD(sn, h), w, sn->onwp[MN_SID_TAIL(sn, t)]);
    });
}

WtU_CHANGESTAT_FN(u__wtsubnets){ 
  GET_AUX_STORAGE(StoreWtSubnets, sn);
  WtSetEdge(MN_IO_TAIL(sn, tail), MN_IO_HEAD(sn, head), weight, sn->onwp[MN_SID_TAIL(sn, tail)]);
}

WtF_CHANGESTAT_FN(f__wtsubnets){
  GET_AUX_STORAGE(StoreWtSubnets, sn);
  for(unsigned int i=1; i<=sn->ns; i++)
    WtNetworkDestroy(sn->onwp[i]);
  sn->onwp++;
  Free(sn->onwp);
}

// MultiNet: Take a weighted networkwise sum of the networks' statistics.

WtI_CHANGESTAT_FN(i_wtMultiNet){
  /*
    iinputs expects:
    1: number of weights (nwts)
    inputs expects:
    nwts*ns: matrix of weights, in network-major order
  */
  
  GET_AUX_STORAGE(StoreWtSubnets, sn);
  unsigned int ns = sn->ns;
  unsigned int nwts = *IINPUT_PARAM;
  double *wts = INPUT_PARAM;
  
  ALLOC_STORAGE(ns, WtModel*, ms);

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
      ms[i-1] = WtModelInitialize(VECTOR_ELT(submodels, submodpos), NULL, sn->onwp[i], FALSE);
      submodpos++;
    }else ms[i-1] = NULL;
  }
  WtDELETE_IF_UNUSED_IN_SUBMODELS(u_func, ms, sn->ns);
  WtDELETE_IF_UNUSED_IN_SUBMODELS(z_func, ms, sn->ns);
}

WtC_CHANGESTAT_FN(c_wtMultiNet){
  GET_AUX_STORAGE(StoreWtSubnets, sn);
  GET_STORAGE(WtModel*, ms);
  unsigned int nwts = *IINPUT_PARAM;
  double *wts = INPUT_PARAM;

  unsigned int i = MN_SID_TAIL(sn, tail);
  WtModel *m = ms[i-1];
  if(m){ // NULL if network has weights 0.
    Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
    WtChangeStats1(st, sh, weight, sn->onwp[i], m, edgeweight);

    wts += (i-1)*nwts; // Position of that network's weight vector.
    for(unsigned int j=0; j<m->n_stats; j++)
      for(unsigned int k=0; k<nwts; k++)
	CHANGE_STAT[j*nwts+k] += m->workspace[j]*wts[k];
  }
}

WtZ_CHANGESTAT_FN(z_wtMultiNet){
  GET_AUX_STORAGE(StoreWtSubnets, sn);
  GET_STORAGE(WtModel*, ms);
  unsigned int nwts = *IINPUT_PARAM;
  double *wts = INPUT_PARAM;

  for(unsigned int i=1; i<=sn->ns; i++){
    WtModel *m = ms[i-1];
    if(m){ // NULL if network has weights 0.
      WtZStats(sn->onwp[i], m, FALSE);
      
      wts += (i-1)*nwts; // Position of that network's weight vector.
      for(unsigned int j=0; j<m->n_stats; j++)
        for(unsigned int k=0; k<nwts; k++)
          CHANGE_STAT[j*nwts+k] += m->workspace[j]*wts[k];
    }
  }
}

WtF_CHANGESTAT_FN(f_wtMultiNet){
  GET_AUX_STORAGE(StoreWtSubnets, sn);
  GET_STORAGE(WtModel*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    if(ms[i-1]) WtModelDestroy(sn->onwp[i], ms[i-1]);
  }
}

// MultiNets: Concatenate the networks' statistics; network statistic counts may be heterogeneous.

WtI_CHANGESTAT_FN(i_wtMultiNets){
  int *iinputs = IINPUT_PARAM; 
  GET_AUX_STORAGE(StoreWtSubnets, sn);
  unsigned int ns = sn->ns;
  unsigned int *pos = (unsigned int *) iinputs;
  ALLOC_STORAGE(ns, WtModel*, ms);

  SEXP submodels = getListElement(mtp->R, "submodels");
  unsigned int submodpos = 0;
  for(unsigned int i=1; i<=sn->ns; i++){
    if(pos[i-1]!=pos[i]){
      ms[i-1] = WtModelInitialize(VECTOR_ELT(submodels, submodpos), NULL, sn->onwp[i], FALSE);
      submodpos++;
    }
  }
  WtDELETE_IF_UNUSED_IN_SUBMODELS(u_func, ms, sn->ns);
  WtDELETE_IF_UNUSED_IN_SUBMODELS(z_func, ms, sn->ns);
}

WtC_CHANGESTAT_FN(c_wtMultiNets){
  unsigned int *pos = (unsigned int *) IINPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreWtSubnets, sn);
  GET_STORAGE(WtModel*, ms);

  unsigned int i = MN_SID_TAIL(sn, tail);
  Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
  if(pos[i-1]!=pos[i]){
    WtModel *m = ms[i-1];
    WtChangeStats1(st, sh, weight, sn->onwp[i], m, edgeweight);
    memcpy(CHANGE_STAT + (unsigned int)(pos[i-1]), m->workspace, m->n_stats*sizeof(double));
  }
}

WtZ_CHANGESTAT_FN(z_wtMultiNets){
  unsigned int *pos = (unsigned int *) IINPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreWtSubnets, sn);
  GET_STORAGE(WtModel*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    if(pos[i-1]!=pos[i]){
      WtModel *m = ms[i-1];
      WtZStats(sn->onwp[i], m, FALSE);
      memcpy(CHANGE_STAT + (unsigned int)(pos[i-1]), m->workspace, m->n_stats*sizeof(double));
    }
  }
}

WtF_CHANGESTAT_FN(f_wtMultiNets){
  unsigned int *pos = (unsigned int *) IINPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreWtSubnets, sn);
  GET_STORAGE(WtModel*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    if(pos[i-1]!=pos[i]){
      WtModelDestroy(sn->onwp[i], ms[i-1]);
    }
  }
}

