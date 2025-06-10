/*  File src/changestats_dgw_sp_ML.c in package ergm.multi, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */
#include "changestats_dgw_sp_ML.h"
#include "changestats.h"

#define all_calcs(term)                         \
  dvec_calc(term ## _ML)                               \
       dist_calc(term ## _ML)                          \
       gw_calc(term ## _ML)

#define all_calcs2(term)                        \
  dvec_calc2(term ## _ML)                              \
       dist_calc2(term ## _ML)                         \
       gw_calc2(term ## _ML)

#define sp_args tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,N_CHANGE_STATS,dvec,CHANGE_STAT

#define dvec_calc(term)                                                 \
  static inline void term ## _calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreStrictDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, Rboolean any_order, int nd, Vertex *dvec, double *cs) { \
    term ## _change({                                                   \
        for(unsigned int j = 0; j < nd; j++){                           \
          Vertex deg = dvec[j];                                         \
          cs[j] += ((L2 + c2path == deg) - (L2 == deg));                \
        }                                                               \
      },{                                                               \
        for(unsigned int j = 0; j < nd; j++){                           \
          Vertex deg = dvec[j];                                         \
          cs[j] += l3c * (L2 == deg);                                   \
        }                                                               \
      });                                                               \
  }

#define dvec_calc2(term)                                                \
  static inline void term ## _calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreStrictDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, Rboolean any_order, int nd, Vertex *dvec, double *cs) { \
    term ## _change({                                                   \
        for(unsigned int j = 0; j < nd; j++){                           \
          Vertex deg = (Vertex)dvec[j];                                 \
          cs[j] += ((L2 + c2path == deg) - (L2 == deg)) * 2;            \
        }                                                               \
      },{                                                               \
        for(unsigned int j = 0; j < nd; j++){                           \
          Vertex deg = (Vertex)dvec[j];                                 \
          cs[j] += l3c * (L2 == deg);                                   \
        }                                                               \
      });                                                               \
  }

#define spd_args tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,N_CHANGE_STATS,CHANGE_STAT

#define dist_calc(term)                                                 \
  static inline void term ## _dist_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreStrictDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, Rboolean any_order, int nd, double *cs) { \
    term ## _change({                                                   \
        int nL2 = L2 + c2path;                                          \
        if(nL2 > nd) cutoff_error(mtp);                                 \
        if(L2) cs[L2-1]--;                                              \
        if(nL2) cs[nL2-1]++;                                            \
      },{                                                               \
        if(L2 > nd) cutoff_error(mtp);                                  \
        if(L2) cs[L2-1] += l3c;                                         \
      });                                                               \
  }

#define dist_calc2(term)                                                \
  static inline void term ## _dist_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreStrictDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, Rboolean any_order, int nd, double *cs) { \
    term ## _change({                                                   \
        int nL2 = L2 + c2path;                                          \
        if(nL2 > nd) cutoff_error(mtp);                                 \
        if(L2) cs[L2-1]-=2;                                             \
        if(nL2) cs[nL2-1]+=2;                                           \
      },{                                                               \
        if(L2 > nd) cutoff_error(mtp);                                  \
        if(L2) cs[L2-1] += l3c;                                         \
      });                                                               \
  }


#define gwsp_args tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,alpha,loneexpa

#define gw_calc(term)                                                   \
  static inline double term ## _gw_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreStrictDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, Rboolean any_order, double alpha, double loneexpa) { \
    double cumchange = 0;                                               \
    term ## _change({                                                   \
        int nL2 = L2 + c2path;                                          \
        cumchange += exp(log1mexp(-loneexpa * nL2)) - exp(log1mexp(-loneexpa * L2)); \
      },{                                                               \
        cumchange += l3c * exp(log1mexp(-loneexpa * L2));               \
      });                                                               \
    return cumchange;                                                   \
  }


#define gw_calc2(term)                                                  \
  static inline double term ## _gw_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreStrictDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, Rboolean any_order, double alpha, double loneexpa) { \
    double cumchange = 0;                                               \
    term ## _change({                                                   \
        int nL2 = L2 + c2path;                                          \
        cumchange += (exp(log1mexp(-loneexpa * nL2)) - exp(log1mexp(-loneexpa * L2))) * 2; \
      },{                                                               \
        cumchange += l3c * exp(log1mexp(-loneexpa * L2));               \
      });                                                               \
    return cumchange;                                                   \
  }


all_calcs(dspUTP)
all_calcs(dspOTP)
all_calcs(dspITP)
all_calcs2(dspOSP)
all_calcs2(dspISP)
/* all_calcs2(dspRTP) */

/*****************
 changestat: d_dsp
*****************/
/*
  Note that d_esp is a meta-function, dispatching actual changescore
  calculation to one of the esp*_calc routines, based on the selected shared
  partner type code.

  Type codes are as follows (where (i,j) is the focal edge):

  UTP - Undirected two-path (undirected graphs only)
  OTP - Outgoing two-path (i->k->j)
  ITP - Incoming two-path (i<-k<-j)
  RTP - Reciprocated two-path (i<->k<->j)
  OSP - Outgoing shared partner (i->k<-j)
  ISP - Incoming shared partner (i<-k->j)

  Only one type may be specified per esp term.  UTP should always be used for undirected graphs; OTP is the traditional directed default.
*/
C_CHANGESTAT_FN(c_ddsp_ML) { 
  /*Set things up*/
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 2);
  StoreLayerLogic *ll3 = NULL;
  StoreStrictDyadMapUInt *spcache = N_AUX>=4 ? AUX_STORAGE_NUM(3) : NULL;
  Rboolean any_order = (Rboolean) IINPUT_PARAM[0];
  unsigned int type=IINPUT_PARAM[1];     /*Get the ESP type code to be used*/
  Vertex *dvec=(Vertex *) IINPUT_PARAM+2;           /*Get the pointer to the ESP stats list*/

  /*Obtain the DSP changescores (by type)*/
  switch(type){
  case ESPUTP: dspUTP_ML_calc(sp_args); break;
  case ESPOTP: dspOTP_ML_calc(sp_args); break;
  case ESPITP: dspITP_ML_calc(sp_args); break;
  /* case ESPRTP: dspRTP_ML_calc(sp_args); break; */
  case ESPOSP: dspOSP_ML_calc(sp_args); break;
  case ESPISP: dspISP_ML_calc(sp_args); break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/  
}


C_CHANGESTAT_FN(c_ddspdist_ML) {
  /*Set things up*/
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 2);
  StoreLayerLogic *ll3 = NULL;
  StoreStrictDyadMapUInt *spcache = N_AUX>=4 ? AUX_STORAGE_NUM(3) : NULL;
  Rboolean any_order = (Rboolean) IINPUT_PARAM[0];
  unsigned int type=IINPUT_PARAM[1];     /*Get the ESP type code to be used*/

  /*Obtain the DSP changescores (by type)*/
  switch(type){
  case ESPUTP: dspUTP_ML_dist_calc(spd_args); break;
  case ESPOTP: dspOTP_ML_dist_calc(spd_args); break;
  case ESPITP: dspITP_ML_dist_calc(spd_args); break;
  /* case ESPRTP: dspRTP_ML_dist_calc(spd_args); break; */
  case ESPOSP: dspOSP_ML_dist_calc(spd_args); break;
  case ESPISP: dspISP_ML_dist_calc(spd_args); break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/
}


/*****************
 changestat: d_gwdsp
*****************/

/*
  Note that d_gwesp is a meta-function for all geometrically weighted ESP stats; the specific type of ESP to be employed is determined by the type argument (INPUT_PARAM[1]).  Type codes are as follows (where (i,j) is the focal edge):

  OTP (0) - Outgoing two-path (i->k->j)
  ITP (1) - Incoming two-path (i<-k<-j)
  RTP (2) - Reciprocated two-path (i<->k<->j)
  OSP (3) - Outgoing shared partner (i->k<-j)
  ISP (4) - Incoming shared partner (i<-k->j)

  Only one type may be specified per esp term.  The default, OTP, retains the original behavior of esp/gwesp.  In the case of undirected graphs, OTP should be used (the others assume a directed network memory structure, and are not safe in the undirected case).
*/
C_CHANGESTAT_FN(c_dgwdsp_ML) {
  /*Set things up*/
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 2);
  StoreLayerLogic *ll3 = NULL;
  StoreStrictDyadMapUInt *spcache = N_AUX>=4 ? AUX_STORAGE_NUM(3) : NULL;
  Rboolean any_order = (Rboolean) IINPUT_PARAM[0];
  double alpha = INPUT_PARAM[0];       /*Get alpha*/
  double loneexpa = log1mexp(alpha);    /*Precompute log(1-exp(-alpha))*/
  int type = IINPUT_PARAM[1];     /*Get the ESP type code to be used*/
  double cumchange = 0;

  /*Obtain the DSP changescores (by type)*/
  switch(type){
    case ESPUTP: cumchange = dspUTP_ML_gw_calc(gwsp_args); break;
    case ESPOTP: cumchange = dspOTP_ML_gw_calc(gwsp_args); break;
    case ESPITP: cumchange = dspITP_ML_gw_calc(gwsp_args); break;
    /* case ESPRTP: dspRTP_ML_gw_calc(gwsp_args); break; */
    case ESPOSP: cumchange = dspOSP_ML_gw_calc(gwsp_args); break;
    case ESPISP: cumchange = dspISP_ML_gw_calc(gwsp_args); break;
  }

  CHANGE_STAT[0] = exp(alpha) * cumchange;
}


all_calcs(espUTP)
all_calcs(espOTP)
all_calcs(espITP)
all_calcs(espOSP)
all_calcs(espISP)
/* all_calcs(espRTP) */


/*****************
 changestat: d_esp
*****************/
/*
  Note that d_esp is a meta-function, dispatching actual changescore
  calculation to one of the esp*_calc routines, based on the selected shared
  partner type code.

  Type codes are as follows (where (i,j) is the focal edge):

  UTP - Undirected two-path (undirected graphs only)
  OTP - Outgoing two-path (i->k->j)
  ITP - Incoming two-path (i<-k<-j)
  RTP - Reciprocated two-path (i<->k<->j)
  OSP - Outgoing shared partner (i->k<-j)
  ISP - Incoming shared partner (i<-k->j)

  Only one type may be specified per esp term.  UTP should always be used for undirected graphs; OTP is the traditional directed default.
*/
C_CHANGESTAT_FN(c_desp_ML) {
  /*Set things up*/
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll3, 3);
  StoreStrictDyadMapUInt *spcache = N_AUX>=5 ? AUX_STORAGE_NUM(4) : NULL;
  Rboolean any_order = (Rboolean) IINPUT_PARAM[0];
  int type = IINPUT_PARAM[1];     /*Get the ESP type code to be used*/
  Vertex *dvec = (Vertex*) IINPUT_PARAM+2;           /*Get the pointer to the ESP stats list*/

  /*Obtain the ESP changescores (by type)*/
  switch(type){
    case ESPUTP: espUTP_ML_calc(sp_args); break;
    case ESPOTP: espOTP_ML_calc(sp_args); break;
    case ESPITP: espITP_ML_calc(sp_args); break;
    /* case ESPRTP: espRTP_ML_calc(sp_args); break; */
    case ESPOSP: espOSP_ML_calc(sp_args); break;
    case ESPISP: espISP_ML_calc(sp_args); break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/
}


C_CHANGESTAT_FN(c_despdist_ML) {
  /*Set things up*/
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll3, 3);
  StoreStrictDyadMapUInt *spcache = N_AUX>=5 ? AUX_STORAGE_NUM(4) : NULL;
  Rboolean any_order = (Rboolean) IINPUT_PARAM[0];
  int type = IINPUT_PARAM[1];     /*Get the ESP type code to be used*/

  /*Obtain the ESP changescores (by type)*/
  switch(type){
  case ESPUTP: espUTP_ML_dist_calc(spd_args); break;
  case ESPOTP: espOTP_ML_dist_calc(spd_args); break;
  case ESPITP: espITP_ML_dist_calc(spd_args); break;
  /* case ESPRTP: espRTP_ML_dist_calc(spd_args); break; */
  case ESPOSP: espOSP_ML_dist_calc(spd_args); break;
  case ESPISP: espISP_ML_dist_calc(spd_args); break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/
}


/*****************
 changestat: d_gwesp
*****************/

/*
  Note that d_gwesp is a meta-function for all geometrically weighted ESP stats; the specific type of ESP to be employed is determined by the type argument (INPUT_PARAM[1]).  Type codes are as follows (where (i,j) is the focal edge):

  OTP (0) - Outgoing two-path (i->k->j)
  ITP (1) - Incoming two-path (i<-k<-j)
  RTP (2) - Reciprocated two-path (i<->k<->j)
  OSP (3) - Outgoing shared partner (i->k<-j)
  ISP (4) - Incoming shared partner (i<-k->j)

  Only one type may be specified per esp term.  The default, OTP, retains the original behavior of esp/gwesp.  In the case of undirected graphs, OTP should be used (the others assume a directed network memory structure, and are not safe in the undirected case).
*/
C_CHANGESTAT_FN(c_dgwesp_ML) { 
  /*Set things up*/
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll3, 3);
  StoreStrictDyadMapUInt *spcache = N_AUX>=5 ? AUX_STORAGE_NUM(4) : NULL;
  Rboolean any_order = (Rboolean) IINPUT_PARAM[0];
  double alpha = INPUT_PARAM[0];       /*Get alpha*/
  double loneexpa = log1mexp(alpha);    /*Precompute (1-exp(-alpha))*/
  int type = IINPUT_PARAM[1];     /*Get the ESP type code to be used*/
  double cumchange = 0;

  /*Obtain the ESP changescores (by type)*/
  switch(type){
    case ESPUTP: cumchange = espUTP_ML_gw_calc(gwsp_args); break;
    case ESPOTP: cumchange = espOTP_ML_gw_calc(gwsp_args); break;
    case ESPITP: cumchange = espITP_ML_gw_calc(gwsp_args); break;
    /* case ESPRTP: cumchange = espRTP_ML_gw_calc(gwsp_args); break; */
    case ESPOSP: cumchange = espOSP_ML_gw_calc(gwsp_args); break;
    case ESPISP: cumchange = espISP_ML_gw_calc(gwsp_args); break;
  }

  CHANGE_STAT[0] = exp(alpha) * cumchange;
}


/*****************
 changestat: d_nsp
*****************/
/*
  Note that d_esp is a meta-function, dispatching actual changescore
  calculation to one of the esp*_calc routines, based on the selected shared
  partner type code.

  Type codes are as follows (where (i,j) is the focal edge):

  UTP - Undirected two-path (undirected graphs only)
  OTP - Outgoing two-path (i->k->j)
  ITP - Incoming two-path (i<-k<-j)
  RTP - Reciprocated two-path (i<->k<->j)
  OSP - Outgoing shared partner (i->k<-j)
  ISP - Incoming shared partner (i<-k->j)

  Only one type may be specified per esp term.  UTP should always be used for undirected graphs; OTP is the traditional directed default.
*/
#define NEGATE_CHANGE_STATS for(unsigned int i = 0; i < N_CHANGE_STATS; i++) CHANGE_STAT[i] *= -1;
C_CHANGESTAT_FN(c_dnsp_ML) {
  /*Set things up*/
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll3, 3);
  StoreStrictDyadMapUInt *spcache = N_AUX>=5 ? AUX_STORAGE_NUM(4) : NULL;
  Rboolean any_order = (Rboolean) IINPUT_PARAM[0];
  int type = IINPUT_PARAM[1];     /*Get the ESP type code to be used*/
  Vertex *dvec = (Vertex*) IINPUT_PARAM+2;           /*Get the pointer to the ESP stats list*/

  /*Obtain the DSP changescores (by type)*/
  switch(type){
  case ESPUTP: 
    espUTP_ML_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspUTP_ML_calc(sp_args);
    break;
  case ESPOTP: 
    espOTP_ML_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspOTP_ML_calc(sp_args);
    break;
  case ESPITP: 
    espITP_ML_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspITP_ML_calc(sp_args);
    break;
  /* case ESPRTP:  */
  /*   espRTP_ML_calc(sp_args); */
  /*   NEGATE_CHANGE_STATS; */
  /*   dspRTP_ML_calc(sp_args); */
  /*   break; */
  case ESPOSP: 
    espOSP_ML_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspOSP_ML_calc(sp_args);
    break;
  case ESPISP: 
    espISP_ML_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspISP_ML_calc(sp_args);
    break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/  
}


C_CHANGESTAT_FN(c_dnspdist_ML) {
  /*Set things up*/
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll3, 3);
  StoreStrictDyadMapUInt *spcache = N_AUX>=5 ? AUX_STORAGE_NUM(4) : NULL;
  Rboolean any_order = (Rboolean) IINPUT_PARAM[0];
  int type = IINPUT_PARAM[1];     /*Get the ESP type code to be used*/

  /*Obtain the DSP changescores (by type)*/
  switch(type){
  case ESPUTP:
    espUTP_ML_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspUTP_ML_dist_calc(spd_args);
    break;
  case ESPOTP:
    espOTP_ML_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspOTP_ML_dist_calc(spd_args);
    break;
  case ESPITP:
    espITP_ML_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspITP_ML_dist_calc(spd_args);
    break;
  /* case ESPRTP:  */
  /*   espRTP_ML_dist_calc(spd_args); */
  /*   NEGATE_CHANGE_STATS; */
  /*   dspRTP_ML_dist_calc(spd_args); */
  /*   break; */
  case ESPOSP:
    espOSP_ML_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspOSP_ML_dist_calc(spd_args);
    break;
  case ESPISP:
    espISP_ML_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspISP_ML_dist_calc(spd_args);
    break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/
}


/*****************
 changestat: d_gwnsp
*****************/

/*
  Note that d_gwesp is a meta-function for all geometrically weighted ESP stats; the specific type of ESP to be employed is determined by the type argument (INPUT_PARAM[1]).  Type codes are as follows (where (i,j) is the focal edge):

  OTP (0) - Outgoing two-path (i->k->j)
  ITP (1) - Incoming two-path (i<-k<-j)
  RTP (2) - Reciprocated two-path (i<->k<->j)
  OSP (3) - Outgoing shared partner (i->k<-j)
  ISP (4) - Incoming shared partner (i<-k->j)

  Only one type may be specified per esp term.  The default, OTP, retains the original behavior of esp/gwesp.  In the case of undirected graphs, OTP should be used (the others assume a directed network memory structure, and are not safe in the undirected case).
*/
C_CHANGESTAT_FN(c_dgwnsp_ML) {
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll3, 3);
  StoreStrictDyadMapUInt *spcache = N_AUX>=5 ? AUX_STORAGE_NUM(4) : NULL;
  Rboolean any_order = (Rboolean) IINPUT_PARAM[0];
  double alpha = INPUT_PARAM[0];       /*Get alpha*/
  double loneexpa = log1mexp(alpha);    /*Precompute (1-exp(-alpha))*/
  int type = IINPUT_PARAM[1];     /*Get the ESP type code to be used*/
  double cumchange = 0;

  /*Obtain the changescores (by type)*/
  switch(type){
  case ESPUTP:
    cumchange = dspUTP_ML_gw_calc(gwsp_args) - espUTP_ML_gw_calc(gwsp_args);
    break;
  case ESPOTP:
    cumchange = dspOTP_ML_gw_calc(gwsp_args) - espOTP_ML_gw_calc(gwsp_args);
    break;
  case ESPITP:
    cumchange = dspITP_ML_gw_calc(gwsp_args) - espITP_ML_gw_calc(gwsp_args);
    break;
  /* case ESPRTP: */
  /*   cumchange = dspRTP_ML_gw_calc(gwsp_args) - espRTP_ML_gw_calc(gwsp_args); */
  /*   break; */
  case ESPOSP:
    cumchange = dspOSP_ML_gw_calc(gwsp_args) - espOSP_ML_gw_calc(gwsp_args);
    break;
  case ESPISP:
    cumchange = dspISP_ML_gw_calc(gwsp_args) - espISP_ML_gw_calc(gwsp_args);
    break;
  }

  CHANGE_STAT[0] = exp(alpha) * cumchange;
}
