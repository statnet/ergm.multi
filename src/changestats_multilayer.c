/*  File src/changestats_multilayer.c in package ergm.multi, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "ergm.multi_changestat_multilayer.h"
#include "ergm_changestat_operator.h"
#include "ergm_changestat.h"
#include "ergm_model.h"
#include "ergm_storage.h"

I_CHANGESTAT_FN(i__layer_net){
  int *iinputs = IINPUT_PARAM;
  ALLOC_AUX_STORAGE(1, StoreLayerLogic, ll);
  ll->nl = *(iinputs++);
  ll->inwp = nwp;

  /* Set up the layer information. */
  ll->lid = (Vertex*) iinputs - 1; // The -1 is because Vertex IDs count from 1.
  iinputs += N_NODES;
  ll->lmap = (Vertex*) iinputs - 1;
  iinputs += N_NODES;
  Vertex lnnodes, lbip;
  if(BIPARTITE){
    lbip = lnnodes = *(iinputs++);
    lnnodes += *(iinputs++);
    iinputs += (ll->nl-1)*2; // There will be a total of nl*2 network sizes.
  }else{
    lbip = 0;
    lnnodes = *(iinputs++);
    iinputs += (ll->nl-1); // There will be a total of nl network sizes.
  }

  if(DIRECTED){
    ll->symm = iinputs - 1; // The -1 is because layer IDs count from 1.
    unsigned int need_symm = FALSE;
    for(unsigned int l=1; l<=ll->nl; l++){
      if(ll->symm[l]){
	need_symm = TRUE;
	break;
      }
    }
    if(!need_symm) ll->symm = NULL;
    
    iinputs += ll->nl;
  }else ll->symm = NULL;

  ll->onwp = NetworkInitialize(NULL, NULL, 0, lnnodes, DIRECTED, lbip);
  
  /* Set up the layer logic. */

  ll->commands = iinputs;
  unsigned int nst = 0;
  int com;
  // com nonnegative means pushing to the stack, either layer or constant.
  for(unsigned int i = 0; (com = ll->commands[i]) != ML_STOP; i++){
    if(com >= 0) nst++;
    if(com == 0) i++; // Next value is a literal.
  }
  ll->stack0 = R_Calloc(nst, int);
  ll->stack1 = R_Calloc(nst, int);

  /* Figure out if this layer needs to calculate reciprocal toggles. */

  ll->need_ht = FALSE;
  if(DIRECTED){
    int com;
    for(unsigned int i = 0; (com = ll->commands[i]) != ML_STOP; i++){
      if(com == 0){ // Next value is a literal.
        i++;
        continue;
      }
      if(com == -21 || // If t() is ever used, or
         (ll->symm && com > 0 && ll->symm[com])){ // any symmetrized layers are referenced,
        ll->need_ht = TRUE; // then toggle (t,h) somewhere may affect dyad (h,t) in this layer.
        break;
      }
    }
  }

  /* Construct the output (logical layer) network: */  

  EXEC_THROUGH_NET_EDGES_PRE(t, h, e, {
      Vertex at[2];
      Vertex ah[2];
      unsigned int nt = ergm_LayerLogic_affects(t, h, ll, LL_ASIS, at, ah);
      for(unsigned int i=0; i<nt; i++){
	ML_SETWT(ll, at[i], ah[i], 1);
      }
    });
}

U_CHANGESTAT_FN(u__layer_net){ 
  GET_AUX_STORAGE(StoreLayerLogic, ll);
  Vertex at[2], ah[2];
  unsigned int nt = ergm_LayerLogic_affects(tail, head, ll, LL_DIFF, at, ah);
  for(unsigned int i=0; i<nt; i++){
    ML_TOGGLE(ll, at[i], ah[i]);
  }
}

F_CHANGESTAT_FN(f__layer_net){ 
  GET_AUX_STORAGE(StoreLayerLogic, ll);
  NetworkDestroy(ll->onwp);
  R_Free(ll->stack0);
  R_Free(ll->stack1);
}

/* I_CHANGESTAT_FN(i__layer_nets){ */
/*   ALLOC_AUX_STORAGE(1, StoreNetsAndLIDAndLMapAndNL, li); */
/*   li->nl = INPUT_PARAM[1]; */
/*   Vertex lnnodes = N_NODES/li->nl, lbip = BIPARTITE/li->nl; */
/*   li->nwp = R_Calloc(li->nl+1, Network); */
/*   li->lid = INPUT_PARAM+2 -1; // The -1 is because Vertex IDs count from 1. */
/*   li->lmap = INPUT_PARAM+2+N_NODES -1; */
/*   for(unsigned int l = 1; l <= li->nl; l++){ */
/*     li->nwp[l] = NetworkInitialize(NULL, NULL, 0, lnnodes, DIRECTED, lbip); */
/*   } */
  
/*   EXEC_THROUGH_NET_EDGES_PRE(t, h, e, { */
/*       ToggleEdge(li->lmap[t], li->lmap[h], li->nwp + (int)li->lid[t]); */
/*     }); */
/* } */

/* U_CHANGESTAT_FN(u__layer_nets){  */
/*   GET_AUX_STORAGE(StoreNetsAndLIDAndLMapAndNL, li); */
/*   ToggleEdge(li->lmap[tail], li->lmap[head], li->nwp + (int)li->lid[tail]); */
/* } */

/* F_CHANGESTAT_FN(f__layer_nets){  */
/*   GET_AUX_STORAGE(StoreNetsAndLIDAndLMapAndNL, li); */
/*   for(unsigned int l = 1; l <= li->nl; l++){ */
/*     NetworkDestroy(li->nwp + l); */
/*   } */
/*   R_Free(li->nwp); */
/* } */

I_CHANGESTAT_FN(i_OnLayer){
  
  unsigned int nml = *IINPUT_PARAM; // Number of layers *in the term*.

  ALLOC_STORAGE(nml, Model*, ms);

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE(ml, StoreLayerLogic, ll);
    ms[ml] = ModelInitialize(getListElement(mtp->R, "submodel"), NULL, ll->onwp, FALSE);
  }
  DELETE_IF_UNUSED_IN_SUBMODELS(u_func, ms, nml);
  DELETE_IF_UNUSED_IN_SUBMODELS(z_func, ms, nml);
}

C_CHANGESTAT_FN(c_OnLayer){
  GET_STORAGE(Model*, ms);
  unsigned int nml = *IINPUT_PARAM;
  double *w = INPUT_PARAM;

  // Find the affected models.
  for(unsigned int ml=0; ml < nml; ml++){
    GET_AUX_STORAGE(ml, StoreLayerLogic, ll);
    Vertex at[2], ah[2];
    unsigned int nt = ergm_LayerLogic_affects(tail, head, ll, LL_DIFF, at, ah);
    if(nt){
      ChangeStats(nt, at, ah, ll->onwp, ms[ml]);
      for(unsigned int i=0; i<N_CHANGE_STATS; i++)
	CHANGE_STAT[i] += ms[ml]->workspace[i] * w[ml];
    }
  }
}

Z_CHANGESTAT_FN(z_OnLayer){
  GET_STORAGE(Model*, ms);
  unsigned int nml = *IINPUT_PARAM;
  double *w = INPUT_PARAM;

  // Find the affected models.
  for(unsigned int ml=0; ml < nml; ml++){
    GET_AUX_STORAGE(ml, StoreLayerLogic, ll);
    ZStats(ll->onwp, ms[ml], FALSE);
    for(unsigned int i=0; i<N_CHANGE_STATS; i++)
      CHANGE_STAT[i] += ms[ml]->workspace[i] * w[ml];
  }
}

F_CHANGESTAT_FN(f_OnLayer){
  GET_STORAGE(Model*, ms);
  unsigned int nml = *IINPUT_PARAM;
  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE(ml, StoreLayerLogic, ll);
    ModelDestroy(ll->onwp, ms[ml]);
  }
}

/* layerCMB: Conway-Maxwell-Binomial for the sum of layer combinations */

C_CHANGESTAT_FN(c_layerCMB){
  unsigned int nml = *IINPUT_PARAM;

  // FIXME: Cache current values, perhaps via a valued auxiliary?

  unsigned int oldct_th=0, newct_th=0,
    oldct_ht=0, newct_ht=0;

  /* Determine whether we need to check the reciprocating dyads. */
  Rboolean need_ht = FALSE;
  for(unsigned int ml=0; ml < nml; ml++){
    GET_AUX_STORAGE(ml, StoreLayerLogic, ll);
    if(ll->need_ht){
      need_ht = TRUE;
      break;
    }
  }

  for(unsigned int ml=0; ml < nml; ml++){
    GET_AUX_STORAGE(ml, StoreLayerLogic, ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    unsigned int v = ergm_LayerLogic2(lt, lh, tail, head, ll, LL_ENCODE);
    if(v&1) oldct_th++; // Pre-toggle edge present.
    if(v&2) newct_th++; // Post-toggle edge present.

    if(need_ht){
      v = ergm_LayerLogic2(lh, lt, tail, head, ll, LL_ENCODE);
      if(v&1) oldct_ht++; // Pre-toggle edge present.
      if(v&2) newct_ht++; // Post-toggle edge present.
    }
  }
  
  CHANGE_STAT[0] =
    +(newct_th!=oldct_th? lgamma1p(newct_th)-lgamma1p(oldct_th) + lgamma1p(nml-newct_th)-lgamma1p(nml-oldct_th) : 0)
    +(newct_ht!=oldct_ht? lgamma1p(newct_ht)-lgamma1p(oldct_ht) + lgamma1p(nml-newct_ht)-lgamma1p(nml-oldct_ht) : 0);    
}

/*****************
 changestat: c_twostarL
*****************/
C_CHANGESTAT_FN(c_twostarL) { 
  unsigned int typeID = IINPUT_PARAM[0];
  unsigned int distinct = IINPUT_PARAM[1];
  
  GET_AUX_STORAGE(0, StoreLayerLogic, ll1);
  GET_AUX_STORAGE(1, StoreLayerLogic, ll2);
  
  Vertex lt = ML_IO_TAIL(ll1, tail), lh = ML_IO_HEAD(ll1, head);

  int
    change1_th = ergm_LayerLogic2(lt, lh, tail, head, ll1, LL_DIFF),
    change1_ht = ergm_LayerLogic2(lh, lt, tail, head, ll1, LL_DIFF),
    change2_th = ergm_LayerLogic2(lt, lh, tail, head, ll2, LL_DIFF),
    change2_ht = ergm_LayerLogic2(lh, lt, tail, head, ll2, LL_DIFF);


  // Need int here since we need signed arithmetic.
  int *od1 = (int*) ML_OUT_DEG(ll1), *od2 = (int*) ML_OUT_DEG(ll2),
    *id1 = (int*) ML_IN_DEG(ll1), *id2 = (int*) ML_IN_DEG(ll2);
  
  switch(typeID){
  case 0: // any / undirected
    if(change1_th || change2_th || change1_ht || change2_ht){ // any counts change.

      // Calculate the change in number of relevant coincident relations.
      int change12 = 0;
      if(distinct){
        if(ML_IS_UNDIRECTED_EDGE(ll1, lt, lh)) change12 += change2_th + change2_ht;
        if(ML_IS_UNDIRECTED_EDGE(ll2, lt, lh)) change12 += change1_th + change1_ht;
        change12 += (change1_th+change1_ht) * (change2_th+change2_ht);
      }

      CHANGE_STAT[0] +=
        + (od1[lt]+id1[lt]+change1_th+change1_ht)*(od2[lt]+id2[lt]+change2_th+change2_ht) - (od1[lt]+id1[lt])*(od2[lt]+id2[lt]) // Change due to lt
        + (od1[lh]+id1[lh]+change1_th+change1_ht)*(od2[lh]+id2[lh]+change2_th+change2_ht) - (od1[lh]+id1[lh])*(od2[lh]+id2[lh]) // Change due to lh

        - change12*2;
    }
    break;
  case 1: // out
    if(change1_th || change2_th){ // lt's outstar counts change.
      
      // Calculate the change in number of relevant coincident relations.
      int change12 = 0;
      if(distinct){
	if(ML_IS_OUTEDGE(ll1, lt, lh)) change12 += change2_th; 
	if(ML_IS_OUTEDGE(ll2, lt, lh)) change12 += change1_th;
	change12 += change1_th * change2_th;
      }
	
      CHANGE_STAT[0] += (od1[lt]+change1_th)*(od2[lt]+change2_th) - od1[lt]*od2[lt] - change12;
    }
    
    if(change1_ht || change2_ht){ // lh's outstar counts change.

      // Calculate the change in number of relevant coincident relations.
      int change12 = 0;
      if(distinct){
	if(ML_IS_OUTEDGE(ll1, lh, lt)) change12 += change2_ht; 
	if(ML_IS_OUTEDGE(ll2, lh, lt)) change12 += change1_ht;
	change12 += change1_ht * change2_ht;
      }

      CHANGE_STAT[0] += (od1[lh]+change1_ht)*(od2[lh]+change2_ht) - od1[lh]*od2[lh] - change12;
    }
    break;
  case 2: // in
    if(change1_th || change2_th){ // lh's instar counts change.
      
      // Calculate the change in number of relevant coincident relations.
      int change12 = 0;
      if(distinct){
	if(ML_IS_OUTEDGE(ll1, lt, lh)) change12 += change2_th; 
	if(ML_IS_OUTEDGE(ll2, lt, lh)) change12 += change1_th;
	change12 += change1_th * change2_th;
      }

      CHANGE_STAT[0] += (id1[lh]+change1_th)*(id2[lh]+change2_th) - id1[lh]*id2[lh] - change12;
    }
    
    if(change1_ht || change2_ht){ // lt's instar counts change.

      // Calculate the change in number of relevant coincident relations.
      int change12 = 0;
      if(distinct){
	if(ML_IS_OUTEDGE(ll1, lh, lt)) change12 += change2_ht; 
	if(ML_IS_OUTEDGE(ll2, lh, lt)) change12 += change1_ht;
	change12 += change1_ht * change2_ht;
      }

      CHANGE_STAT[0] += (id1[lt]+change1_ht)*(id2[lt]+change2_ht) - id1[lt]*id2[lt] - change12;
    }
    break;
  case 3: // path
    if(change1_ht || change2_th){ // two-path count through lt changes.

      // Calculate the change in number of relevant coincident relations.
      int change12 = 0;
      if(distinct){
	if(ML_IS_OUTEDGE(ll1, lh, lt)) change12 += change2_th; 
	if(ML_IS_OUTEDGE(ll2, lt, lh)) change12 += change1_ht;
	change12 += change1_th * change2_ht;
      }

      CHANGE_STAT[0] += (id1[lt]+change1_ht)*(od2[lt]+change2_th) - id1[lt]*od2[lt] - change12;
    }
    
    if(change1_th || change2_ht){ // two-path count through lh changes.

      // Calculate the change in number of relevant coincident relations.
      int change12 = 0;
      if(distinct){
	if(ML_IS_INEDGE(ll1, lh, lt)) change12 += change2_ht; 
	if(ML_IS_OUTEDGE(ll2, lh, lt)) change12 += change1_th;
	change12 += change1_th * change2_ht;
      }

      CHANGE_STAT[0] += (id1[lh]+change1_th)*(od2[lh]+change2_ht) - id1[lh]*od2[lh] - change12;
    }
    break;
  }
}

/*****************
 changestat: d_mutual_ML

 (1,1) -> anything = -1
 anything -> (1,1) = +1
*****************/
C_CHANGESTAT_FN(c_mutual_ML){
  GET_AUX_STORAGE(0, StoreLayerLogic, ll1);
  GET_AUX_STORAGE(1, StoreLayerLogic, ll2);
  
  double matchval;
  int j, ninputs, noattr;

  ninputs = N_INPUT_PARAMS - N_NODES;
  noattr = !N_INPUT_PARAMS;

  /* *** don't forget tail -> head */
  Vertex lt = ML_IO_TAIL(ll1, tail), lh = ML_IO_HEAD(ll1, head);
  int l1th = ergm_LayerLogic2(lt, lh, tail, head, ll1, LL_ENCODE);
  int l1ht = ergm_LayerLogic2(lh, lt, tail, head, ll1, LL_ENCODE);
  int l2th = ergm_LayerLogic2(lt, lh, tail, head, ll2, LL_ENCODE);
  int l2ht = ergm_LayerLogic2(lh, lt, tail, head, ll2, LL_ENCODE);

  int change =
    +((l1th&2)&&(l2ht&2))-((l1th&1)&&(l2ht&1)) // t-l1->h and h->l2->t
    +((l2th&2)&&(l1ht&2))-((l2th&1)&&(l1ht&1)) // t-l2->h and h->l1->t
    ;
  
  if(change) { /* otherwise, no change occurs */
      if (noattr) { /* "plain vanilla" mutual, without node attributes */
        CHANGE_STAT[0] += change;
      } else { /* Only consider mutuals where node attributes match */
        matchval = INPUT_PARAM[tail+ninputs-1];
        if (matchval == INPUT_PARAM[head+ninputs-1]) { /* We have a match! */
          if (ninputs==0) {/* diff=F in network statistic specification */
            CHANGE_STAT[0] += change;
          } else { /* diff=T */
            for (j=0; j<ninputs; j++) {
              if (matchval == INPUT_PARAM[j]) 
                CHANGE_STAT[j] += change;
            }
          }
        }
      }
    }
}

/*****************
 changestat: d_mutual_by_attr_ML
*****************/
C_CHANGESTAT_FN(c_mutual_by_attr_ML) { 
  GET_AUX_STORAGE(0, StoreLayerLogic, ll1);
  GET_AUX_STORAGE(1, StoreLayerLogic, ll2);

  int j, ninputs;

  ninputs = N_INPUT_PARAMS - N_NODES;

  /* *** don't forget tail -> head */
  Vertex lt = ML_IO_TAIL(ll1, tail), lh = ML_IO_HEAD(ll1, head);
  int l1th = ergm_LayerLogic2(lt, lh, tail, head, ll1, LL_ENCODE);
  int l1ht = ergm_LayerLogic2(lh, lt, tail, head, ll1, LL_ENCODE);
  int l2th = ergm_LayerLogic2(lt, lh, tail, head, ll2, LL_ENCODE);
  int l2ht = ergm_LayerLogic2(lh, lt, tail, head, ll2, LL_ENCODE);

  int change =
    +((l1th&2)&&(l2ht&2))-((l1th&1)&&(l2ht&1)) // t-l1->h and h->l2->t
    +((l2th&2)&&(l1ht&2))-((l2th&1)&&(l1ht&1)) // t-l2->h and h->l1->t
    ;

    if (change) { /* otherwise, no change occurs */
      for (j=0; j<ninputs; j++) {
        if (INPUT_PARAM[tail+ninputs-1] == INPUT_PARAM[j]){CHANGE_STAT[j] += change;}
        if (INPUT_PARAM[head+ninputs-1] == INPUT_PARAM[j]){CHANGE_STAT[j] += change;}
      }
    }
}

typedef enum {LAYER_DIST_HAMMING = 0, LAYER_DIST_ENTRAINMENT = 1} LayerDistType;

/*****************
 changestat: c_pairwisedistL
*****************/
C_CHANGESTAT_FN(c_pairwisedistL){
  /* Obtain the list of layers (to be compared) that are potentially
     affected by the toggle. */
  LayerDistType type = *IINPUT_PARAM;
  Vertex *rdeps = ((Vertex *) IINPUT_PARAM) + 1;
  StoreLayerLogic *ll0 = AUX_STORAGE;
  Vertex tl = ML_LID_TAIL(ll0, tail), nl = N_AUX;
  Vertex *affected = rdeps + rdeps[tl - 1],
    n_affected = rdeps[tl] - rdeps[tl - 1];

  /* If none affected, no change. */
  if(n_affected == 0) return;

  /* Find out if we need reciprocal ties. */
  Rboolean need_ht = FALSE;
  for(Vertex i = 0; i < n_affected; i++) {
    Vertex l = affected[i];
    GET_AUX_STORAGE(l, StoreLayerLogic, ll);
    if(ll->need_ht){
      need_ht = TRUE;
      break;
    }
  }

  Rboolean th_before[nl], th_after[nl],
    ht_before[need_ht ? nl : 1], ht_after[need_ht ? nl : 1];

  /* Evaluate old states. */
  for(Vertex l = 0; l < nl; l++) {
    GET_AUX_STORAGE(l, StoreLayerLogic, ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    th_before[l] = ML_GETWT(ll, lt, lh);
    if(need_ht) ht_before[l] = ML_GETWT(ll, lh, lt);
  }

  memcpy(th_after, th_before, nl * sizeof(Rboolean));
  if(need_ht) memcpy(ht_after, ht_before, nl * sizeof(Rboolean));

  /* Update new states where applicable. */
  for(Vertex i = 0; i < n_affected; i++) {
    Vertex l = affected[i];
    GET_AUX_STORAGE(l, StoreLayerLogic, ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    th_after[l] = ergm_LayerLogic2(lt, lh, tail, head, ll, LL_POST);
    if(need_ht) ht_after[l] = ergm_LayerLogic2(lh, lt, tail, head, ll, LL_POST);
  }

  /* Calculate the changes in Hamming distances. */
  for(Vertex i = 0; i < n_affected; i++) {
    Vertex l1 = affected[i];
    if(th_after[l1] != th_before[l1])
      for(Vertex l2 = 0; l2 < nl; l2++)
        if(l1 != l2)
          switch(type){
          case LAYER_DIST_HAMMING:
            CHANGE_STAT[0] += (th_after[l1] != th_after[l2]) - (th_before[l1] != th_before[l2]);
            break;
          case LAYER_DIST_ENTRAINMENT:
            CHANGE_STAT[0] += (th_after[l1] && th_after[l2]) - (th_before[l1] && th_before[l2]);
            break;
          }

    if(need_ht && ht_after[l1] != ht_before[l1])
      for(Vertex l2 = 0; l2 < nl; l2++)
        if(l1 != l2)
          switch(type){
          case LAYER_DIST_HAMMING:
            CHANGE_STAT[0] += (ht_after[l1] != ht_after[l2]) - (ht_before[l1] != ht_before[l2]);
            break;
          case LAYER_DIST_ENTRAINMENT:
            CHANGE_STAT[0] += (ht_after[l1] && ht_after[l2]) - (ht_before[l1] && ht_before[l2]);
            break;
          default: error("This should not happen.");
          }
  }
}
