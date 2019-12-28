/*  File src/MHProposals_block.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#include "ergm_MHproposal.h"
#include "ergm_edgelist.h"
#include "ergm_changestat.h"
#include "ergm_MHblockdiag.h"
#include "ergm_MHstorage.h"

/*********************
 void MH_blockdiag

 Block-diagonal sampling
*********************/
MH_I_FN(Mi_blockdiag){
  ALLOC_STORAGE(1, MH_BlockDiagSampInfo, b);
  *b = unpack_BlockDiagSampInfo(getListElement(MHp->R,"BDI"), BIPARTITE, DIRECTED);
  MHp->ntoggles=1;
}

MH_P_FN(Mp_blockdiag){
  GET_STORAGE(MH_BlockDiagSampInfo, b);

  BD_LOOP({
      GetRandDyadBlockDiag(Mtail, Mhead, b);
    });
}

/********************
   void MH_blockTNT

   Block-diagonal TNT sampling
***********************/
MH_I_FN(Mi_blockdiagTNT){
  ALLOC_STORAGE(1, MH_BlockDiagSampInfo, b);
  *b = unpack_BlockDiagSampInfo(getListElement(MHp->R,"BDI"), BIPARTITE, DIRECTED);
  MHp->ntoggles=1;
}

MH_P_FN(Mp_blockdiagTNT){
  GET_STORAGE(MH_BlockDiagSampInfo, b);

  const double P = 0.5, Q = 1-P;
  double DP = P*b->ndyads, DO = DP/Q;

  Edge nedges=EDGECOUNT(nwp);
  
  double logratio=0; 

  BD_LOOP({
      if (unif_rand() < P && nedges > 0) { /* Select a tie at random */
	// Note that, by construction, this tie will be within a block.
	GetRandEdge(Mtail, Mhead, nwp);
	logratio = TNT_LR_E(nedges, Q, DP, DO);
      }else{ /* Select a dyad at random within a block */
	GetRandDyadBlockDiag(Mtail, Mhead, b);
	
	if(IS_OUTEDGE(Mtail[0],Mhead[0])!=0){
          logratio = TNT_LR_DE(nedges, Q, DP, DO);
	}else{
          logratio = TNT_LR_DN(nedges, Q, DP, DO);
	}
      }
    });
  MHp->logratio += logratio;
}
