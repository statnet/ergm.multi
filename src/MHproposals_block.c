/*  File src/MHproposals_block.c in package ergm.multi, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#include "ergm_MHproposal.h"
#include "ergm_MHproposal_bd.h"
#include "ergm_edgelist.h"
#include "ergm_changestat.h"
#include "ergm.multi_MHblockdiag.h"
#include "ergm_MHstorage.h"

typedef struct{
  MH_BlockDiagSampInfo b;
  DegreeBound *bd;
} StoreBlockDiagSampInfoAndDegreeBound;

/*********************
 void MH_blockdiag

 Block-diagonal sampling
*********************/
MH_I_FN(Mi_blockdiag){
  ALLOC_STORAGE(1, StoreBlockDiagSampInfoAndDegreeBound, storage);
  storage->b = unpack_BlockDiagSampInfo(getListElement(MHp->R,"BDI"), BIPARTITE, DIRECTED);
  storage->bd = DegreeBoundInitializeR(MHp->R, nwp);
  MHp->ntoggles=1;
}

MH_P_FN(Mp_blockdiag){
  GET_STORAGE(StoreBlockDiagSampInfoAndDegreeBound, storage);

  BD_LOOP(storage->bd, {
      GetRandDyadBlockDiag(Mtail, Mhead, &storage->b);
    });
}

MH_F_FN(Mf_blockdiag) {
  GET_STORAGE(StoreBlockDiagSampInfoAndDegreeBound, storage);

  DegreeBoundDestroy(storage->bd);
}


/********************
   void MH_blockTNT

   Block-diagonal TNT sampling
***********************/
MH_I_FN(Mi_blockdiagTNT){
  ALLOC_STORAGE(1, StoreBlockDiagSampInfoAndDegreeBound, storage);
  storage->b = unpack_BlockDiagSampInfo(getListElement(MHp->R,"BDI"), BIPARTITE, DIRECTED);
  storage->bd = DegreeBoundInitializeR(MHp->R, nwp);
  MHp->ntoggles=1;
}

MH_P_FN(Mp_blockdiagTNT){
  GET_STORAGE(StoreBlockDiagSampInfoAndDegreeBound, storage);

  const double P = 0.5, Q = 1-P;
  double DP = P*storage->b.ndyads, DO = DP/Q;

  Edge nedges=EDGECOUNT(nwp);
  
  double logratio=0; 

  BD_LOOP(storage->bd, {
      if (unif_rand() < P && nedges > 0) { /* Select a tie at random */
	// Note that, by construction, this tie will be within a block.
	GetRandEdge(Mtail, Mhead, nwp);
	logratio = TNT_LR_E(nedges, Q, DP, DO);
      }else{ /* Select a dyad at random within a block */
	GetRandDyadBlockDiag(Mtail, Mhead, &storage->b);
	
	if(IS_OUTEDGE(Mtail[0],Mhead[0])!=0){
          logratio = TNT_LR_DE(nedges, Q, DP, DO);
	}else{
          logratio = TNT_LR_DN(nedges, Q, DP, DO);
	}
      }
    });
  MHp->logratio += logratio;
}

MH_F_FN(Mf_blockdiagTNT) {
  GET_STORAGE(StoreBlockDiagSampInfoAndDegreeBound, storage);

  DegreeBoundDestroy(storage->bd);
}
