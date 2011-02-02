#include "wtMHproposals.h"
/* Shorthand. */
#define Mhead (MHp->togglehead)
#define Mtail (MHp->toggletail)
#define Mweight (MHp->toggleweight)

/*********************
 void MH_Poisson

 Default MH algorithm for Poisson-reference ERGM
*********************/
void MH_BipartitePoisson(WtMHproposal *MHp, WtNetwork *nwp)  {  
  Vertex head, tail;
  double oldwt;
  int fvalid, trytoggle;
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=1;
    return;
  }
  MHp->ratio = 1.0;
  
  Mhead[0] = 1 + unif_rand() * nwp->bipartite;
  Mtail[0] = 1 + nwp->bipartite + unif_rand() * (nwp->nnodes - nwp->bipartite);

  oldwt = WtGetEdge(Mhead[0],Mtail[0],nwp);

  const double fudge = 0.5; // Mostly comes in when proposing from 0.

  do{
    Mweight[0] = rpois(oldwt + fudge);    
  }while(Mweight[0]==oldwt);
    
  MHp->ratio *= exp((1 + log(Mweight[0]+fudge))*oldwt - (1 + log(oldwt+fudge))*Mweight[0]) * (1-dpois(oldwt,oldwt+fudge,0))/(1-dpois(Mweight[0],Mweight[0]+fudge,0));
}

/*********************
 void MH_CompleteOrderingBipartite

 Default MH algorithm for ERGM over complete orderings
*********************/
void MH_CompleteOrderingBipartite(WtMHproposal *MHp, WtNetwork *nwp)  {  
  Vertex head, tail1, tail2;
  double oldwt;
  int fvalid, trytoggle;
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=2;
    return;
  }
    
  
  Mhead[0] = Mhead[1] = 1 + unif_rand() * nwp->bipartite;
  Mtail[0] = 1 + nwp->bipartite + unif_rand() * (nwp->nnodes - nwp->bipartite);
  Mtail[1] = 1 + Mtail[0] + unif_rand() * (nwp->nnodes - Mtail[0]);
  
  Mweight[1] = WtGetEdge(Mhead[0],Mtail[0],nwp);
  Mweight[0] = WtGetEdge(Mhead[1],Mtail[1],nwp);
}
