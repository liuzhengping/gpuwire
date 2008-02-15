#ifndef _DELTA_KERNEL_H_
#define _DELTA_KERNEL_H_

#include <stdio.h>
#define DELTA 1.0
#define INF 1e20
#define BUCKETSIZE 4096*8*2
#define NUMBUCKETS 512

#define DOWN  0
#define UP    1
#define RIGHT 2
#define LEFT  3


//#define EMULATION

#ifdef EMULATION
#define DEBUG(x...) printf(x)
#else
#define DEBUG(x...) 
#endif

texture<float, 2, cudaReadModeElementType> mytex0;
texture<float, 2, cudaReadModeElementType> mytex1;
texture<float, 2, cudaReadModeElementType> mytex2;
texture<float, 2, cudaReadModeElementType> mytex3;

__device__ void
demptyKernel(){
}

__global__ void
labelKernel (int i, int* B,int* BCount,int* BPos, int tw,int th, int* RLoc,int* R, float* dR,int* dSourceR, float* d,int* vBucketMap){
//ver se não precisa colocar 4 BiCount em algum lugar...

  //todo: try to increase speed using shared memory for RLoc... think more about it (maybe RLoc is too big for shared memory)
  const unsigned int tid = threadIdx.x;
  const unsigned int num_threads = blockDim.x;

  int BiCount = BPos[i];
  int node, row, col,index;
  float cost,f1,f2,fmin;

  //cleaning R
  //4 times because each node can be reached from up, down, left and right directions (and more 4 times because for each node 4 more are open)
//  DEBUG("BiCount %d lastpos %d\n",BiCount,4*(BiCount-1)+3);
  for(int k=0;  (num_threads*k + tid) < 16*BiCount;k++){
    R[(num_threads*k + tid)  ]=-1;
    dR[(num_threads*k + tid)  ]=INF;
    dSourceR[(num_threads*k + tid)  ]=-1;
  }


  __syncthreads();
  const int dx[4]={0,0,1,-1};
  const int dy[4]={1,-1,0,0};
  

  for(int k=0; num_threads*k + tid < BiCount;k++){
    node = B  [ BUCKETSIZE*i + num_threads*k + tid];
    if(node!=-1){

      DEBUG("(tid %d) node %d(from B[%d] pos %d)\n",tid,node,i,num_threads*k + tid);
      for(int j=0;j<4;j++){
        switch(j){
	  case 0:
           cost = tex2D(mytex0,node%tw,node/tw);
           break;
	  case 1:
           cost = tex2D(mytex1,node%tw,node/tw);
           break;
	  case 2:
           cost = tex2D(mytex2,node%tw,node/tw);
           break;
	  case 3:
           cost = tex2D(mytex3,node%tw,node/tw);
           break;
        }
        row = node/tw + dy[j];
        col = node%tw + dx[j];

        if( (row>=0) && (row < th) && (col >= 0) && (col < tw) ){
          RLoc [ row*tw + col ] = 4*(num_threads*k + tid)+j;
          DEBUG("(tid %d)Connecting node %d to be processed by %d\n",tid,row*tw+col,RLoc [ row*tw + col ]);
        }

      }

    }
  }

  __syncthreads();

  //copy Edges to R
  for(int k=0; num_threads*k + tid < BiCount;k++){
    node = B  [ BUCKETSIZE*i + num_threads*k + tid];
    if(node!=-1){
      for(int j=0;j<4;j++){
        switch(j){
  	case 0:
         cost = tex2D(mytex0,node%tw,node/tw);
         break;
	case 1:
         cost = tex2D(mytex1,node%tw,node/tw);
         break;
	case 2:
         cost = tex2D(mytex2,node%tw,node/tw);
         break;
	case 3:
         cost = tex2D(mytex3,node%tw,node/tw);
         break;
      }



        row = node/tw + dy[j];
        col = node%tw + dx[j];
              DEBUG("Pre-candidate in R %d (d=%f)\n",row*tw + col,d[node]+cost,dR [4*RLoc[row*tw + col]+j]);
        if( (row>=0) && (row < th) && (col >= 0) && (col < tw) ){ 
            DEBUG("Candidate in R %d (d=%f)\n",row*tw + col,d[node]+cost,dR [4*RLoc[row*tw + col]+j]);
          if((cost<=DELTA)&&( d[node]+cost < d[row*tw + col])){
              R [4*RLoc[row*tw + col]+j] = row*tw + col;
              dR[4*RLoc[row*tw + col]+j] = d[node]+cost;
	      dSourceR[4*RLoc[row*tw + col]+j]=node;
              DEBUG("New node in R %d (d=%f,e=%f) in pos %d\n",R [4*RLoc[row*tw + col]+j],dR [4*RLoc[row*tw + col]+j],cost,4*RLoc[row*tw + col]+j);
              vBucketMap[node]=-1;
          }
        }
      }
    }
  }
  __syncthreads();

  //gathering data to find the minimum cost way to get to node n
  //TODO: OPTIMIZE IN SUCH A WAY IT WON'T BE NEEDED TO GO THROUGH THE 4 EDGES, since they store the same value
  int smin ;
  float dists[4];
  for(int k=0;  (num_threads*k + tid) < 4*BiCount;k++){
    index    = 4*(num_threads*k + tid);
    dists[0] = dR[ index  ];
    dists[1] = dR[ index+1];
    dists[2] = dR[ index+2];
    dists[3] = dR[ index+3];

    //finds the node with minimum distance, so that the path can be stored
    if( dists[0] < dists[1]){
	if( dists[2] < dists[3]){
	  if( dists[0] < dists[2]){
	    smin = dSourceR[index];
	  }
          else{
            smin = dSourceR[index+2];
	  }
	}
        else{
	  if( dists[0] < dists[3]){
	    smin = dSourceR[index];
	  }
	  else{
	    smin = dSourceR[index+3];
	  }
	}
    }
    else{
	if( dists[2] < dists[3]){
	  if( dists[1] < dists[2]){
	    smin = dSourceR[index+1];
	  }
          else{
            smin = dSourceR[index+2];
	  }
	}
        else{
	  if( dists[1] < dists[3]){
	    smin = dSourceR[index+1];
	  }
	  else{
	    smin = dSourceR[index+3];
	  }
	}

    }    

    f1 = fminf( dists[0] , dists[1] );
    f2 = fminf( dists[2] , dists[3] );

    fmin = fminf(f1,f2);
    dR[index  ]=fmin;
    dR[index+1]=fmin;
    dR[index+2]=fmin;
    dR[index+3]=fmin;

    DEBUG("Smin %d\n",smin);
    dSourceR[index  ]=smin;
    dSourceR[index+1]=smin;
    dSourceR[index+2]=smin;
    dSourceR[index+3]=smin;

    
  }
  __syncthreads();

}

//Pensar se o fato de S ter duplicatas pode impactar em algo
__global__ void 
copyB2SKernel(int i, int* B, int* BCount,int* BPos, int* S, int* SCount){
  //TODO: optimize this code
  //there's an optimized way of doing this, which is by only
  //storing SCount = Scount+ BCount, as output
  //and controlling with local variables thread positions

  const unsigned int tid = threadIdx.x;
  const unsigned int num_threads = blockDim.x;

  int pos;

  int BiCount = BCount[i];
  for(int k=0; num_threads*k + tid < BiCount;k++){
    if(B[i*BUCKETSIZE+num_threads*k+tid]!=-1){
      pos = atomicAdd(&SCount[0],1);
      S[pos] = B[i*BUCKETSIZE+num_threads*k+tid];
    }
  }
  __syncthreads();
  BCount[i]=0;
  BPos[i]=0;
  __syncthreads();
//  DEBUG("(tid %d) SCount %d\n",tid,SCount[0]);


}


//Parallel relax edges
__global__ void
relaxKernelPath( int RCount, int* B,int* BCount,int* BPos, int* RLoc,int* R,float*  dR,float* d, int* vBucketLoc, int* vBucketMap, float* deb, int* dSourceR, int* dSource){

  const unsigned int tid = threadIdx.x;
  const unsigned int num_threads = blockDim.x;
  int v,bn,bn_old, index;
  float x;


  //remove node from old bucket



  for(int k=0; num_threads*k + tid < RCount;k++){
    index = num_threads*k + tid;

    if(R[index]!=-1){

      x = dR[index];      
      v = R[index];

      if(x<d[v]){

        bn_old = vBucketMap[v];
        if (bn_old != -1) {
          int oldIndex = bn_old*BUCKETSIZE+vBucketLoc[v];
//	  
          B[oldIndex] = -1;//GN;
//          int oldc = atomicSub(&BCount[bn_old],1);
	  atomicSub(&BCount[bn_old],1);

//        printf("Removing %d from %d(%d)\n",v,bn_old,oldc);


        }
      }

    }    

  }




  __syncthreads();



  for(int k=0; num_threads*k + tid < RCount;k++){

    if(R[num_threads*k + tid]!=-1){
//deb[0]= (float) (BUCKETSIZE);
      x = dR[num_threads*k + tid];
      v = R[num_threads*k + tid];
      if(x < d[v]){

	      bn = (int) (dR[num_threads*k + tid]/DELTA);

	      atomicAdd(&BCount[bn],1);
	      int pos = atomicAdd(&BPos[bn],1);
	      DEBUG("Pos %d BCount[%d] %d node %d (x=%f)\n",pos,bn,BPos[bn],v,x);
	      
	      B[bn*BUCKETSIZE+pos] = v;
	      d[v] = x;
	      dSource[v] = dSourceR[num_threads*k+tid];
	      
	      vBucketLoc[v] = pos;
	      vBucketMap[v] = bn;
	      RLoc[v]=-1;
      }

    }
  }


  __syncthreads();

  
}



//Parallel relax edges
__global__ void
relaxKernel( int RCount, int* B,int* BCount,int* BPos, int* RLoc,int* R,float*  dR,float* d, int* vBucketLoc, int* vBucketMap, float* deb){

  const unsigned int tid = threadIdx.x;
  const unsigned int num_threads = blockDim.x;
  int v,bn,bn_old, index;
  float x;
//  int myAdd=0;
  deb[20]= (float)RCount;

  for(int k=0; num_threads*k + tid < RCount;k++){
    index = num_threads*k + tid;
    deb[2*index]= R[index];
    deb[2*index+1]= dR[index];
  }

//  DEBUG("relaxing RCount %d\n",RCount);

  //remove node from old bucket
//  RCount = RCount /4;



  for(int k=0; num_threads*k + tid < RCount;k++){
    index = num_threads*k + tid;
  //  deb[index]= R[index];

    if(R[index]!=-1){

      x = dR[index];      
      v = R[index];

      if(x<d[v]){

        bn_old = vBucketMap[v];
        if (bn_old != -1) {
          int oldIndex = bn_old*BUCKETSIZE+vBucketLoc[v];
//	  
          B[oldIndex] = -1;//GN;
//          int oldc = atomicSub(&BCount[bn_old],1);
	  atomicSub(&BCount[bn_old],1);

//        printf("Removing %d from %d(%d)\n",v,bn_old,oldc);


        }
      }

    }    

  }




  __syncthreads();






  for(int k=0; num_threads*k + tid < RCount;k++){

    if(R[num_threads*k + tid]!=-1){
//deb[0]= (float) (BUCKETSIZE);
      x = dR[num_threads*k + tid];
      v = R[num_threads*k + tid];
      if(x < d[v]){



      bn = (int) (dR[num_threads*k + tid]/DELTA); 
      
   //   printf("Bn %d\n",bn);

      atomicAdd(&BCount[bn],1);
      int pos = atomicAdd(&BPos[bn],1);
      DEBUG("Pos %d BCount[%d] %d node %d (x=%f)\n",pos,bn,BPos[bn],v,x);

      
      B[bn*BUCKETSIZE+pos] = v;
      d[v] = x;

      
      vBucketLoc[v] = pos;  
      vBucketMap[v] = bn;
      RLoc[v]=-1;


  //only debug info
//    for(int i=0;i<BPos[bn];i++){
//      DEBUG("B(%d)=%d ",i,B[bn*BUCKETSIZE+i]);
//    }
//    DEBUG("\n");



      }

    }
  }


  __syncthreads();

  
}


__global__ void
labelHeavyKernel (int i, int* B,int SCount, int tw,int th, int* RLoc,int* R, float* dR, float* d,int* vBucketMap){
  //todo: try to increase speed using shared memory for RLoc... think more about it (maybe RLoc is too big for shared memory)
  const unsigned int tid = threadIdx.x;
  const unsigned int num_threads = blockDim.x;

  int BiCount = SCount;
  int node, row, col,index;
  float cost,f1,f2,fmin;

  //cleaning R
  //4 times because each node can be reached from up, down, left and right directions (and more 4 times because for each node 4 more are open)
  //  DEBUG("BiCount %d lastpos %d\n",BiCount,4*(BiCount-1)+3);
  for(int k=0;  (num_threads*k + tid) < 16*BiCount;k++){
    index = 4*(num_threads*k + tid);
    R[index  ]=-1;
    R[index+1]=-1;
    R[index+2]=-1;
    R[index+3]=-1;
  
    dR[index  ]=INF;
    dR[index+1]=INF;
    dR[index+2]=INF;
    dR[index+3]=INF;
  }

  __syncthreads();
  const int dx[4]={0,0,1,-1};
  const int dy[4]={1,-1,0,0};  

  for(int k=0; num_threads*k + tid < BiCount;k++){
    node = B  [ BUCKETSIZE*i + num_threads*k + tid];
    if(node!=-1){

      DEBUG("(tid %d) node %d(from B[%d] pos %d)\n",tid,node,i,num_threads*k + tid);
      for(int j=0;j<4;j++){
        //TRY TO OPTIMIZE THIS PART REMOVING THE IF AND SETTING IMPOSSIBLE EDGE VALUES TO INFINITY
        switch(j){
	  case 0:
           cost = tex2D(mytex0,node%tw,node/tw);
           break;
	  case 1:
           cost = tex2D(mytex1,node%tw,node/tw);
           break;
	  case 2:
           cost = tex2D(mytex2,node%tw,node/tw);
           break;
	  case 3:
           cost = tex2D(mytex3,node%tw,node/tw);
           break;
        }
//      cost = tex2D(mytex[j],node%tw,node/tw);//todo: change texture
        row = node/tw + dy[j];
        col = node%tw + dx[j];

        if( (row>=0) && (row < th) && (col >= 0) && (col < tw) ){
          RLoc [ row*tw + col ] = 4*(num_threads*k + tid)+j;
          DEBUG("(tid %d)Connecting node %d to be processed by %d\n",tid,row*tw+col,RLoc [ row*tw + col ]);
        }

      }

    }
  }

  __syncthreads();

  //copy Edges to R
  for(int k=0; num_threads*k + tid < BiCount;k++){
    node = B  [ BUCKETSIZE*i + num_threads*k + tid];
    if(node!=-1){
      for(int j=0;j<4;j++){
        switch(j){
  	case 0:
         cost = tex2D(mytex0,node%tw,node/tw);
         break;
	case 1:
         cost = tex2D(mytex1,node%tw,node/tw);
         break;
	case 2:
         cost = tex2D(mytex2,node%tw,node/tw);
         break;
	case 3:
         cost = tex2D(mytex3,node%tw,node/tw);
         break;
      }



        row = node/tw + dy[j];
        col = node%tw + dx[j];
              DEBUG("Pre-candidate in R %d (d=%f)\n",row*tw + col,d[node]+cost,dR [4*RLoc[row*tw + col]+j]);
        if( (row>=0) && (row < th) && (col >= 0) && (col < tw) ){ 
            DEBUG("Candidate in R %d (d=%f)\n",row*tw + col,d[node]+cost,dR [4*RLoc[row*tw + col]+j]);
          if((cost>DELTA)&&( d[node]+cost < d[row*tw + col])){
              R [4*RLoc[row*tw + col]+j] = row*tw + col;
              dR[4*RLoc[row*tw + col]+j] = d[node]+cost;
              DEBUG("New node in R %d (d=%f,e=%f) in pos %d\n",R [4*RLoc[row*tw + col]+j],dR [4*RLoc[row*tw + col]+j],cost,4*RLoc[row*tw + col]+j);
              vBucketMap[node]=-1;
          }
        }
      }
    }
  }
  __syncthreads();

  //gathering data to find the minimum cost way to get to node n
  //TODO: OPTIMIZE IN SUCH A WAY IT WON'T BE NEEDED TO GO THROUGH THE 4 EDGES, since they store the same value
  for(int k=0;  (num_threads*k + tid) < 16*BiCount;k++){
    f1 = fminf( dR[4*(num_threads*k + tid)  ], dR[4*(num_threads*k + tid)+1] );
    f2 = fminf( dR[4*(num_threads*k + tid)+2], dR[4*(num_threads*k + tid)+3] );
    fmin = fminf(f1,f2);
    dR[4*(num_threads*k + tid)  ]=fmin;
    dR[4*(num_threads*k + tid)+1]=fmin;
    dR[4*(num_threads*k + tid)+2]=fmin;
    dR[4*(num_threads*k + tid)+3]=fmin;
  }
  __syncthreads();
  
}



__global__ void
emptyKernel(){
}
// B is the bucket i vector
// RLoc[n] stores the position of node n in R (so that if more than one attempt to update
// the distance to node n is made at the same time, it can be shifted to 0,1,2 or 3 in the position of R)


#endif // #ifndef _MEMORY_KERNEL_H_
