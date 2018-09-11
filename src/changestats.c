/*  File src/changestats.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_changestat_multilayer.h"
#include "ergm_storage.h"
#include "ergm_dyad_hashmap.h"

/*****************
 changestat: d_b1degree_ML_sum
*****************/
C_CHANGESTAT_FN(c_b1degree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double *degs = inputs;

  /* *** don't forget tail -> head */

  unsigned int b1deg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    b1deg += od[lt];
  }
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex deg = (Vertex)degs[j];
    CHANGE_STAT[j] += (b1deg + degchange == deg) - (b1deg == deg);
  }
}

/*****************
 changestat: d_b1degree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_b1degree_by_attr_ML_sum) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);

  /* *** don't forget tail -> head */    

  unsigned int b1deg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *od=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    b1deg += od[lt];
  }   

  int tailattr = inputs[2*N_CHANGE_STATS + head - 1]; 
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex d = (Vertex)inputs[2*j];
    int testattr = inputs[2*j + 1]; 
    if (tailattr == testattr)  /* we have head attr match */
      CHANGE_STAT[j] += (b1deg + degchange == d) - (b1deg == d);
  }
}

/*****************
 changestat: d_b2degree_ML_sum
*****************/
C_CHANGESTAT_FN(c_b2degree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double *degs = inputs;

  /* *** don't forget tail -> head */

  unsigned int b2deg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    b2deg += id[lh];
  }
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex deg = (Vertex)degs[j];
    CHANGE_STAT[j] += (b2deg + degchange == deg) - (b2deg == deg);
  }
}

/*****************
 changestat: d_b2degree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_b2degree_by_attr_ML_sum) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);

  /* *** don't forget tail -> head */    

  unsigned int b2deg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    b2deg += id[lh];
  }   

  int headattr = inputs[2*N_CHANGE_STATS + head - 1]; 
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex d = (Vertex)inputs[2*j];
    int testattr = inputs[2*j + 1]; 
    if (headattr == testattr)  /* we have head attr match */
      CHANGE_STAT[j] += (b2deg + degchange == d) - (b2deg == d);
  }
}

/*****************
 changestat: d_degree_ML_sum
*****************/
C_CHANGESTAT_FN(c_degree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double *degs = inputs;

  /* *** don't forget tail -> head */

  unsigned int taildeg = 0, headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll), *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    taildeg += od[lt] + id[lt];
    headdeg += od[lh] + id[lh];
  }
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex deg = (Vertex)degs[j];
    CHANGE_STAT[j] += (taildeg + degchange == deg) - (taildeg == deg);
    CHANGE_STAT[j] += (headdeg + degchange == deg) - (headdeg == deg);
  }
}

/*****************
 changestat: d_degree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_degree_by_attr_ML_sum) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);

  /* *** don't forget tail -> head */    

  unsigned int taildeg = 0, headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll), *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    taildeg += od[lt] + id[lt];
    headdeg += od[lh] + id[lh];
  }   

  int tailattr = inputs[2*N_CHANGE_STATS + tail - 1]; 
  int headattr = inputs[2*N_CHANGE_STATS + head - 1]; 
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex d = (Vertex)inputs[2*j];
    int testattr = inputs[2*j + 1]; 
    if (tailattr == testattr)  /* we have tail attr match */
      CHANGE_STAT[j] += (taildeg + degchange == d) - (taildeg == d);
    if (headattr == testattr)  /* we have head attr match */
      CHANGE_STAT[j] += (headdeg + degchange == d) - (headdeg == d);
  }
}

/*****************
 changestat: d_gwb1degree_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwb1degree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double oneexpd = 1.0-exp(-decay);

  unsigned int b1deg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    b1deg += od[lt];
  }

    
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] =
    exp(decay) * ((1-pow(oneexpd,b1deg+degchange)) - (1-pow(oneexpd,b1deg)));
}

/*****************
 changestat: d_gwb1degree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwb1degree_by_attr_ML_sum) { 
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 200?)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double *attrs = inputs-1;
  double oneexpd = 1.0-exp(-decay);

  unsigned int b1deg = 0, headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    b1deg += od[lt];
  }

    
  /* *** don't forget tail -> head */    
  CHANGE_STAT[(int)attrs[tail]] =
    exp(decay) * ((1-pow(oneexpd,b1deg+degchange)) - (1-pow(oneexpd,b1deg)));
}

/*****************
 changestat: d_gwdegree_ML
*****************/
C_CHANGESTAT_FN(c_gwdegree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double oneexpd = 1.0-exp(-decay);

  unsigned int taildeg = 0, headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll), *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    taildeg += od[lt] + id[lt];
    headdeg += od[lh] + id[lh];
  }

  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] =
    exp(decay) * (
		  ((1-pow(oneexpd,taildeg+degchange)) - (1-pow(oneexpd,taildeg))) +
		  ((1-pow(oneexpd,headdeg+degchange)) - (1-pow(oneexpd,headdeg)))
		  );
}

/*****************
 changestat: d_gwdegree_by_attr
*****************/
C_CHANGESTAT_FN(c_gwdegree_by_attr_ML_sum) { 
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 200?)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double *attrs = inputs-1;
  double oneexpd = 1.0-exp(-decay);

  unsigned int taildeg = 0, headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll), *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    taildeg += od[lt] + id[lt];
    headdeg += od[lh] + id[lh];
  }

    
  /* *** don't forget tail -> head */    
  CHANGE_STAT[(int)attrs[tail]] =
    exp(decay) * ((1-pow(oneexpd,taildeg+degchange)) - (1-pow(oneexpd,taildeg)));
  
  CHANGE_STAT[(int)attrs[head]] =
    exp(decay) * ((1-pow(oneexpd,headdeg+degchange)) - (1-pow(oneexpd,headdeg)));
}

/*****************
 changestat: d_gwb2degree_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwb2degree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double oneexpd = 1.0-exp(-decay);

  unsigned int b2deg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    b2deg += id[lh];
  }

    
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] =
        exp(decay) * ((1-pow(oneexpd,b2deg+degchange)) - (1-pow(oneexpd,b2deg)));
}

/*****************
 changestat: d_gwb2degree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwb2degree_by_attr_ML_sum) { 
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 200?)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double *attrs = inputs-1;
  double oneexpd = 1.0-exp(-decay);

  unsigned int b2deg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    b2deg += id[lh];
  }

    
  /* *** don't forget tail -> head */    
  CHANGE_STAT[(int)attrs[head]] =
    exp(decay) * ((1-pow(oneexpd,b2deg+degchange)) - (1-pow(oneexpd,b2deg)));
}

/*****************
 changestat: d_gwidegree_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwidegree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double oneexpd = 1.0-exp(-decay);

  unsigned int headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    headdeg += id[lh];
  }

    
  /* *** don't forget tail -> head */
  CHANGE_STAT[0] =
    exp(decay) * ((1-pow(oneexpd,headdeg+degchange)) - (1-pow(oneexpd,headdeg)));      
}

/*****************
 changestat: d_gwidegree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwidegree_by_attr_ML_sum) { 
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 200?)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double *attrs = inputs-1;
  double oneexpd = 1.0-exp(-decay);

  unsigned int headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    headdeg += id[lh];
  }

    
  /* *** don't forget tail -> head */    
  CHANGE_STAT[(int)attrs[head]] =
    exp(decay) * ((1-pow(oneexpd,headdeg+degchange)) - (1-pow(oneexpd,headdeg)));      
}

/*****************
 changestat: d_gwodegree_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwodegree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double oneexpd = 1.0-exp(-decay);

  unsigned int taildeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    taildeg += od[lt];
  }
    
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] =
    exp(decay) * ((1-pow(oneexpd,taildeg+degchange)) - (1-pow(oneexpd,taildeg)));      
}

/*****************
 changestat: d_gwodegree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwodegree_by_attr_ML_sum) { 
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 200?)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double *attrs = inputs-1;
  double oneexpd = 1.0-exp(-decay);

  unsigned int taildeg = 0, headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    taildeg += od[lt];
  }

    
  /* *** don't forget tail -> head */    
  CHANGE_STAT[(int)attrs[tail]] =
    exp(decay) * ((1-pow(oneexpd,taildeg+degchange)) - (1-pow(oneexpd,taildeg)));
}

/*****************
 changestat: d_idegree_ML_sum
*****************/
C_CHANGESTAT_FN(c_idegree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double *degs = inputs;

  /* *** don't forget tail -> head */

  unsigned int headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    headdeg += id[lh];
  }
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex deg = (Vertex)degs[j];
    CHANGE_STAT[j] += (headdeg + degchange == deg) - (headdeg == deg);
  }
}

/*****************
 changestat: d_degree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_idegree_by_attr_ML_sum) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);

  /* *** don't forget tail -> head */    

  unsigned int headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    headdeg += id[lh];
  }   

  int headattr = inputs[2*N_CHANGE_STATS + head - 1]; 
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex d = (Vertex)inputs[2*j];
    int testattr = inputs[2*j + 1]; 
    if (headattr == testattr)  /* we have head attr match */
      CHANGE_STAT[j] += (headdeg + degchange == d) - (headdeg == d);
  }
}

/*****************
 changestat: d_mutual_ML

 (1,1) -> anything = -1
 anything -> (1,1) = +1
*****************/
C_CHANGESTAT_FN(c_mutual_ML){
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 1);
  
  double matchval;
  int j, ninputs, noattr;

  ninputs = N_INPUT_PARAMS - N_NODES - 2;
  noattr = (N_INPUT_PARAMS == 2);

  /* *** don't forget tail -> head */
  Vertex lt = ML_IO_TAIL(ll1, tail), lh = ML_IO_HEAD(ll1, head);
  int l1th = ergm_LayerLogic2(lt, lh, tail, head, ll1, 2);
  int l1ht = ergm_LayerLogic2(lh, lt, tail, head, ll1, 2);
  int l2th = ergm_LayerLogic2(lt, lh, tail, head, ll2, 2);
  int l2ht = ergm_LayerLogic2(lh, lt, tail, head, ll2, 2);

  int change =
    +((l1th&2)&&(l2ht&2))-((l1th&1)&&(l2ht&1)) // t-l1->h and h->l2->t
    +((l2th&2)&&(l1ht&2))-((l2th&1)&&(l1ht&1)) // t-l2->h and h->l1->t
    ;
  
  if(change) { /* otherwise, no change occurs */
      if (noattr) { /* "plain vanilla" mutual, without node attributes */
        CHANGE_STAT[0] += change;
      } else { /* Only consider mutuals where node attributes match */
        matchval = INPUT_PARAM[tail+ninputs-1+2];
        if (matchval == INPUT_PARAM[head+ninputs-1+2]) { /* We have a match! */
          if (ninputs==0) {/* diff=F in network statistic specification */
            CHANGE_STAT[0] += change;
          } else { /* diff=T */
            for (j=0; j<ninputs; j++) {
              if (matchval == INPUT_PARAM[j+2]) 
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
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 1);

  int j, ninputs;

  ninputs = N_INPUT_PARAMS - N_NODES - 2;

  /* *** don't forget tail -> head */
  Vertex lt = ML_IO_TAIL(ll1, tail), lh = ML_IO_HEAD(ll1, head);
  int l1th = ergm_LayerLogic2(lt, lh, tail, head, ll1, 2);
  int l1ht = ergm_LayerLogic2(lh, lt, tail, head, ll1, 2);
  int l2th = ergm_LayerLogic2(lt, lh, tail, head, ll2, 2);
  int l2ht = ergm_LayerLogic2(lh, lt, tail, head, ll2, 2);

  int change =
    +((l1th&2)&&(l2ht&2))-((l1th&1)&&(l2ht&1)) // t-l1->h and h->l2->t
    +((l2th&2)&&(l1ht&2))-((l2th&1)&&(l1ht&1)) // t-l2->h and h->l1->t
    ;

    if (change) { /* otherwise, no change occurs */
      for (j=0; j<ninputs; j++) {
        if (INPUT_PARAM[tail+ninputs-1+2] == INPUT_PARAM[j+2]){CHANGE_STAT[j] += change;}
        if (INPUT_PARAM[head+ninputs-1+2] == INPUT_PARAM[j+2]){CHANGE_STAT[j] += change;}
      }
    }
}

/*****************
 changestat: d_odegree_ML_sum
*****************/
C_CHANGESTAT_FN(c_odegree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double *degs = inputs;

  /* *** don't forget tail -> head */

  unsigned int taildeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    taildeg += od[lt];
  }
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex deg = (Vertex)degs[j];
    CHANGE_STAT[j] += (taildeg + degchange == deg) - (taildeg == deg);
  }
}

/*****************
 changestat: d_odegree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_odegree_by_attr_ML_sum) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);

  /* *** don't forget tail -> head */    

  unsigned int taildeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    taildeg += od[lt];
  }   

  int tailattr = inputs[2*N_CHANGE_STATS + tail - 1]; 
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex d = (Vertex)inputs[2*j];
    int testattr = inputs[2*j + 1]; 
    if (tailattr == testattr)  /* we have tail attr match */
      CHANGE_STAT[j] += (taildeg + degchange == d) - (taildeg == d);
  }
}
