#include "ergm.multi_changestat_multinet.h"
#include "ergm_changestat_operator.h"
#include "ergm_changestat.h"
#include "ergm_model.h"
#include "ergm_storage.h"


I_CHANGESTAT_FN(i__subnets){
  int *iinputs = IINPUT_PARAM;
  ALLOC_AUX_STORAGE(1, StoreSubnets, sn);
  sn->ns = *(iinputs++);
  sn->inwp = nwp;
  sn->onwp = Calloc(sn->ns, Network *);
  sn->onwp--; // The -- is because Network IDs count from 1.

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

    sn->onwp[i] = NetworkInitialize(NULL, NULL, 0, lnnodes, DIRECTED, lbip, 0, 0, NULL);
  }
  
  EXEC_THROUGH_NET_EDGES_PRE(t, h, e, {
      ToggleKnownEdge(MN_IO_TAIL(sn, t), MN_IO_HEAD(sn, h), sn->onwp[MN_SID_TAIL(sn, t)], FALSE);
    });
}

U_CHANGESTAT_FN(u__subnets){ 
  GET_AUX_STORAGE(StoreSubnets, sn);
  ToggleKnownEdge(MN_IO_TAIL(sn, tail), MN_IO_HEAD(sn, head),sn->onwp[MN_SID_TAIL(sn, tail)], edgeflag);
}

F_CHANGESTAT_FN(f__subnets){
  GET_AUX_STORAGE(StoreSubnets, sn);
  for(unsigned int i=1; i<=sn->ns; i++)
    NetworkDestroy(sn->onwp[i]);
  sn->onwp++;
  Free(sn->onwp);
}

// MultiNet: Take a weighted networkwise sum of the networks' statistics.

I_CHANGESTAT_FN(i_MultiNet){
  /*
    iinputs expects:
    1: number of weights (nwts)
    inputs expects:
    nwts*ns: matrix of weights, in network-major order
  */
  
  GET_AUX_STORAGE(StoreSubnets, sn);
  unsigned int ns = sn->ns;
  unsigned int nwts = *IINPUT_PARAM;
  double *wts = INPUT_PARAM;
  
  ALLOC_STORAGE(ns, Model*, ms);

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
      ms[i-1] = ModelInitialize(VECTOR_ELT(submodels, submodpos), NULL, sn->onwp[i], FALSE);
      submodpos++;
    }else ms[i-1] = NULL;
  }
  DELETE_IF_UNUSED_IN_SUBMODELS(u_func, ms, sn->ns);
  DELETE_IF_UNUSED_IN_SUBMODELS(z_func, ms, sn->ns);
}

C_CHANGESTAT_FN(c_MultiNet){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);
  unsigned int nwts = *IINPUT_PARAM;
  double *wts = INPUT_PARAM;

  unsigned int i = MN_SID_TAIL(sn, tail);
  Model *m = ms[i-1];
  if(m){ // NULL if network has weights 0.
    Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
    ChangeStats1(st, sh, sn->onwp[i], m, edgeflag);

    wts += (i-1)*nwts; // Position of that network's weight vector.
    for(unsigned int j=0; j<m->n_stats; j++)
      for(unsigned int k=0; k<nwts; k++)
	CHANGE_STAT[j*nwts+k] += m->workspace[j]*wts[k];
  }
}

U_CHANGESTAT_FN(u_MultiNet){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);

  unsigned int i = MN_SID_TAIL(sn, tail);
  Model *m = ms[i-1];
  if(m){ // NULL if network has weights 0.
    Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
    UPDATE_STORAGE(st, sh, sn->onwp[i], m, NULL, edgeflag);
  }
}

F_CHANGESTAT_FN(f_MultiNet){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    if(ms[i-1]) ModelDestroy(sn->onwp[i], ms[i-1]);
  }
}

Z_CHANGESTAT_FN(z_MultiNet){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);
  unsigned int nwts = *IINPUT_PARAM;
  double *wts = INPUT_PARAM;

  for(unsigned int i=1; i<=sn->ns; i++){
    Model *m = ms[i-1];
    if(m){ // NULL if network has weights 0.
      ZStats(sn->onwp[i], m, FALSE);

      wts += (i-1)*nwts; // Position of that network's weight vector.
      for(unsigned int j=0; j<m->n_stats; j++)
        for(unsigned int k=0; k<nwts; k++)
          CHANGE_STAT[j*nwts+k] += m->workspace[j]*wts[k];
    }
  }
}

// MultiNets: Concatenate the networks' statistics; network statistic counts may be heterogeneous.

I_CHANGESTAT_FN(i_MultiNets){
  int *iinputs = IINPUT_PARAM; 
  GET_AUX_STORAGE(StoreSubnets, sn);
  unsigned int ns = sn->ns;
  unsigned int *pos = (unsigned int *) iinputs;
  ALLOC_STORAGE(ns, Model*, ms);

  SEXP submodels = getListElement(mtp->R, "submodels");
  unsigned int submodpos = 0;
  for(unsigned int i=1; i<=sn->ns; i++){
    if(pos[i-1]!=pos[i]){
      ms[i-1] = ModelInitialize(VECTOR_ELT(submodels, submodpos), NULL, sn->onwp[i], FALSE);
      submodpos++;
    }
  }
  DELETE_IF_UNUSED_IN_SUBMODELS(u_func, ms, sn->ns);
  DELETE_IF_UNUSED_IN_SUBMODELS(z_func, ms, sn->ns);
}

C_CHANGESTAT_FN(c_MultiNets){
  unsigned int *pos = (unsigned int *) IINPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);

  unsigned int i = MN_SID_TAIL(sn, tail);
  Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
  if(pos[i-1]!=pos[i]){
    Model *m = ms[i-1];
    ChangeStats1(st, sh, sn->onwp[i], m, edgeflag);
    memcpy(CHANGE_STAT + (unsigned int)(pos[i-1]), m->workspace, m->n_stats*sizeof(double));
  }
}

U_CHANGESTAT_FN(u_MultiNets){
  unsigned int *pos = (unsigned int *) IINPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);

  unsigned int i = MN_SID_TAIL(sn, tail);
  if(pos[i-1]!=pos[i]){
    Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
    UPDATE_STORAGE(st, sh, sn->onwp[i], ms[i-1], NULL, edgeflag);
  }
}

Z_CHANGESTAT_FN(z_MultiNets){
  unsigned int *pos = (unsigned int *) IINPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    if(pos[i-1]!=pos[i]){
      Model *m = ms[i-1];
      ZStats(sn->onwp[i], m, FALSE);
      memcpy(CHANGE_STAT + (unsigned int)(pos[i-1]), m->workspace, m->n_stats*sizeof(double));
    }
  }
}

F_CHANGESTAT_FN(f_MultiNets){
  unsigned int *pos = (unsigned int *) IINPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    if(pos[i-1]!=pos[i]){
      ModelDestroy(sn->onwp[i], ms[i-1]);
    }
  }
}

// ByNetDStats

I_CHANGESTAT_FN(i_ByNetDStats){
  Model *m = STORAGE = ModelInitialize(getListElement(mtp->R, "submodel"), NULL, nwp, FALSE);
  DELETE_IF_UNUSED_IN_SUBMODEL(u_func, m);
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);
}

C_CHANGESTAT_FN(c_ByNetDStats){
  unsigned int *pos = (unsigned int *) IINPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model, m);

  unsigned int i = MN_SID_TAIL(sn, tail);
  if(pos[i-1]!=pos[i]){
    ChangeStats1(tail, head, nwp, m, edgeflag);
    memcpy(CHANGE_STAT + (unsigned int)pos[i], m->workspace, m->n_stats*sizeof(double));
  }
}

Z_CHANGESTAT_FN(z_ByNetDStats){
  unsigned int *pos = (unsigned int *) IINPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model, m);

  for(unsigned int i=1; i<=sn->ns; i++)
    if(pos[i-1]!=pos[i]){
      ZStats(nwp, m, FALSE);
      memcpy(CHANGE_STAT + (unsigned int)pos[i], m->workspace, m->n_stats*sizeof(double));
    }
}

U_CHANGESTAT_FN(u_ByNetDStats){
  GET_STORAGE(Model, m);
  UPDATE_STORAGE(tail, head, nwp, m, NULL, edgeflag);
}

F_CHANGESTAT_FN(f_ByNetDStats){
  GET_STORAGE(Model, m);
  ModelDestroy(nwp, m);
  STORAGE = NULL;
}
