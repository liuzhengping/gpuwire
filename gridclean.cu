//compile with
//nvcc delta.cu -I /home/baggio/NVIDIA_CUDA_SDK/common/inc/ -L /home/baggio/NVIDIA_CUDA_SDK/lib/ -lcuda -lcudart -lcutil -lGL -lGLU
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

// includes, project
#include <cutil.h>
#include <cuda_gl_interop.h>
//#define TIMER
#define FINETUNING
// includes, kernels
#include <delta_kernel.cu>
#define NUMTHREADS 128





float* pixels = NULL;
float* dEdges[4];
float maxDistance = 0.0f;
int* hSource = NULL;

int GN = 262144;
int maxRCount = 0;

cudaArray* array[4];


struct edge{
  int dNode[4];
  float weight[4];
};

edge* nodes;





void
runTest( int argc, char** argv,int iw, int ih, int startNode);

void 
runTestGraph ( int argc, char** argv, int iw, int ih, int startNode, float* gradient);


void printPath(int* source, int destination){
	printf("%d <- %d\n",destination,source[destination]);
	if(source[destination]==destination) return;
	printPath( source, source[destination]);
	
}

void printPathSource(int destination, int* path){
	path[0] = hSource[destination];
	printf("%d <- %d\n",destination,hSource[destination]);
	if(hSource[destination]==destination) return;
	printPathSource( hSource[destination],path+1);
}

void loadTexture(int iw, int ih, float* data, cudaArray* cArray, texture<float, 2, cudaReadModeElementType>* myTex){
    cudaChannelFormatDesc desc = cudaCreateChannelDesc<float>();
    CUDA_SAFE_CALL(cudaMallocArray(&cArray, &desc, iw, ih));
    CUDA_SAFE_CALL(cudaMemcpyToArray(cArray, 0, 0, data, sizeof(float)*iw*ih, cudaMemcpyHostToDevice));
    // Bind the array to the texture
    cudaBindTextureToArray( *myTex, cArray, desc);

}

void loadGraphEdges(int iw, int ih,char* myFile){
  int n;
  FILE* in = fopen(myFile,"r");
  fscanf(in,"%d\n",&n);
  
  for(int i=0;i<4;i++){
    dEdges[i] = (float*) malloc(iw*ih*sizeof(float));
    for(int j=0;j<iw*ih;j++){
      dEdges[i][j]=INF;
    }
  }

  while(1){
    int source, dest;
    double eWeight;
    fscanf(in,"%d",&source);
    //printf("Reading %d\n",source);
    if(source==-1) break;
    fscanf(in,"%d%lf\n",&dest,&eWeight);


    //links from 512.graph are either right or down
    if(dest == source + 1){ //it's a RIGHT link
      dEdges[RIGHT][source] = eWeight;
      dEdges[LEFT] [dest  ] = eWeight;
    }
    else if(dest == source + iw){ //it's a DOWN link
      dEdges[DOWN][source] = eWeight;
      dEdges[UP]  [dest  ] = eWeight;
    }        

  }

    loadTexture(iw,ih,dEdges[0],array[0],&mytex0);
    loadTexture(iw,ih,dEdges[1],array[1],&mytex1);
    loadTexture(iw,ih,dEdges[2],array[2],&mytex2);
    loadTexture(iw,ih,dEdges[3],array[3],&mytex3);

  
}

void printDistances(float* hDist,int iw, int ih){
   for(int i=0;i<iw*ih;i++){
      if(i%iw==0) printf("%d ",i/iw);
      if(hDist[i]<INF){
        printf("%5.1f ",hDist[i]);
	if(hDist[i]>maxDistance) maxDistance=hDist[i];
      }
      else
        printf("INFINI ",hDist[i]);
      if(i%iw==iw-1) printf("\n");
    }      
    return;
}



void
runTest( int argc, char** argv, int iw, int ih, int startNode) 
{



    //initialize the device
    cudaSetDevice(0);

    GN = iw*ih;


    unsigned int num_threads = NUMTHREADS;

    // setup execution parameters
    dim3  grid( 1, 1, 1);
    dim3  threads( num_threads, 1, 1);

    nodes = (edge*) malloc(GN*sizeof(edge));   

    //loadGraphEdges(iw,ih,argv[1]);


    float* dDist;
    cudaMalloc( (void**) &dDist, GN*sizeof(float));
    float* hDist = (float*) malloc(GN*sizeof(float));
    
    for(int i=0;i<GN;i++){
      hDist[i]=INF;
    }
    int* dSource;
    cudaMalloc( (void**) &dSource, GN*sizeof(int));
    hSource = (int*) malloc(GN*sizeof(float));
    for(int i=0;i<GN;i++){
      hSource[i]=-1;
    }



    
    int* dBucketMap;
    cudaMalloc( (void**) &dBucketMap, GN*sizeof(int));
    int* hBucketMap = (int*) malloc(GN*sizeof(int));
    
     
    for(int i=0;i<GN;i++){
      hBucketMap[i]=-1;
    }

    cudaMemcpy( dBucketMap, hBucketMap, GN*sizeof(int), cudaMemcpyHostToDevice);

    int* dBucketPos;
    cudaMalloc( (void**) &dBucketPos, GN*sizeof(int));
    int* hBucketPos = (int*) malloc(GN*sizeof(int));

    int* dB;
    cudaMalloc( (void**) &dB, BUCKETSIZE*NUMBUCKETS*sizeof(int));
    int* hB = (int*) malloc(BUCKETSIZE*NUMBUCKETS*sizeof(int));
  

    printf("Bucket size %d\n",BUCKETSIZE);

    int* hBi = hB;
//    int hBiCount[1];
//    hBiCount[0] = 1;
    for(int i=0;i<BUCKETSIZE*NUMBUCKETS;i++)
      hBi[i]=-1;
    
    //start node

    hBi[0] = startNode;
    hDist[startNode] = 0.0f;

    //set startNode source node in the path as the node itself, so that the getPath function works
    cudaMemcpy( dSource+startNode, &startNode, 1*sizeof(int), cudaMemcpyHostToDevice);
    
    



    // copy host memory to device
    cudaMemcpy( dDist, hDist, GN*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy( dB, hBi, BUCKETSIZE*NUMBUCKETS*sizeof(int), cudaMemcpyHostToDevice);


    int* BCount;
    cudaMalloc( (void**) &BCount, NUMBUCKETS * sizeof(int));
    int* hBCount;
    hBCount = (int*) malloc(NUMBUCKETS*sizeof(int));
    for(int i=0;i<NUMBUCKETS;i++)
      hBCount[i]=0;

    hBCount[0]=1;

    // copy host memory to device
    cudaMemcpy( BCount, hBCount, NUMBUCKETS*sizeof(int), cudaMemcpyHostToDevice);


    int* dBPos;
    cudaMalloc( (void**) &dBPos, NUMBUCKETS * sizeof(int));
    int* hBPos;
    hBPos = (int*) malloc(NUMBUCKETS*sizeof(int));
    for(int i=0;i<NUMBUCKETS;i++)
      hBPos[i]=0;
    hBPos[0]=1;
    cudaMemcpy( dBPos, hBPos, NUMBUCKETS*sizeof(int), cudaMemcpyHostToDevice);




    // allocate device memory for result
    int * dRLoc;
    cudaMalloc( (void**) &dRLoc, GN*sizeof(int));
    int * dR;
    cudaMalloc( (void**) &dR, 16*BUCKETSIZE*sizeof(int)); //R is 4 times bigger than RLoc, because UP,DOWN,LEFT and RIGHT will each have a pos in R

    int* hR = (int*) malloc(16*BUCKETSIZE*sizeof(int));
    int* hRLoc = (int*) malloc(GN*sizeof(int));
    float* hDistR= (float*) malloc(16*BUCKETSIZE*sizeof(float));


    float* dDistR;
    cudaMalloc( (void**) &dDistR  , 16*BUCKETSIZE*sizeof(float));
    int* dSourceR;
    cudaMalloc( (void**) &dSourceR, 16*BUCKETSIZE*sizeof(int  ));


    int * dS;
    cudaMalloc( (void**) &dS, 32*BUCKETSIZE*sizeof(int));
    int * hS= (int*) malloc(32*BUCKETSIZE*sizeof(int));
    int * dSCount;
    cudaMalloc( (void**) &dSCount, 1*sizeof(int));
    int * hSCount= (int*) malloc(1*sizeof(int));

    printf("Starting timer\n");

  unsigned int nvtimer = 0;
    cutCreateTimer( &nvtimer);
    cutStartTimer( nvtimer);


  unsigned int laptimer = 1;
    cutCreateTimer( &laptimer);
  float ktime;
   
  unsigned int labelTimer =2;
    cutCreateTimer( &labelTimer);
  int labelCounter=0;

  unsigned int copyTimer =3;
    cutCreateTimer( &copyTimer);
  int copyCounter=0;
    
  unsigned int relaxTimer =4;
    cutCreateTimer( &relaxTimer);
  int relaxCounter=0;


  unsigned int emptyTimer =5;
    cutCreateTimer( &emptyTimer);
  int emptyCounter=0;

  

       float* lido;
    cudaMalloc( (void**) &lido, GN*sizeof(float));
     float* Hlido = (float*) malloc(GN*sizeof(float)) ;
    

      int* RCount;
      RCount = (int*) malloc(1*sizeof(int));
  
    for(int i=0;i<470;i++){
      printf("Bucket %d\n",i);
#ifdef TIMER
      printf("%d\n",i);
#endif

      cudaMemcpy( RCount,  &dBPos[i],  1*sizeof(int),       cudaMemcpyDeviceToHost) ;    
   

/*	printf("Before, Rcount[0] was %d\n",RCount[0]);
      for(int j=0;j<30;j++){
	printf("Bcount[%d]=%d\n",j,debugVector[j]);
      } */


      //S <- EMPTY
      hSCount[0]=0;
      cudaMemcpy( dSCount, hSCount,1*sizeof(int), cudaMemcpyHostToDevice);
    
//  printf("Before RCount %d(i=%d)\n",RCount[0],i);
  //    printf("Still safe%d\n",i);

      //While B[i] != EMPTY
      int sameCount = 0;
      cutStartTimer( laptimer);
      while(RCount[0]!=0){
        if(RCount[0]>maxRCount) maxRCount = RCount[0];
	printf("RCount %d\n",RCount[0]);
//      if(RCount[0]==0) i++;
        sameCount++;

#ifdef TIMER
        cutStartTimer( laptimer);
#endif

#ifdef FINETUNING
        cutStartTimer( labelTimer);
#endif
	
        labelKernel    <<<grid, threads >>> ( i, dB, BCount,dBPos, iw,ih, dRLoc, dR, dDistR, dSourceR, dDist,dBucketMap);
        cudaThreadSynchronize();
    


#ifdef FINETUNING
        cutStopTimer( labelTimer);
	labelCounter++;
#endif

        cudaMemcpy( &hBPos[i],  &dBPos[i],  1*sizeof(int),       cudaMemcpyDeviceToHost) ;


#ifdef TIMER
        ktime = cutGetTimerValue( laptimer );
        printf("Label         kernel %f\n",ktime);	
        cutStartTimer( laptimer);
#endif

#ifdef FINETUNING
        cutStartTimer( copyTimer);
#endif



        copyB2SKernel  <<<grid, threads >>> ( i, dB, BCount,dBPos, dS, dSCount);
        cudaThreadSynchronize();

#ifdef FINETUNING
        cutStopTimer( copyTimer);
	copyCounter++;
#endif



#ifdef TIMER
        ktime = cutGetTimerValue( laptimer );
        printf("CopyB2S       kernel %f\n",ktime);	
#endif

#ifdef TIMER
        cutStartTimer( laptimer);
#endif


#ifdef FINETUNING
        cutStartTimer( emptyTimer);
#endif


//        emptyKernel    <<<grid, threads >>> ();
//        cudaThreadSynchronize();

#ifdef FINETUNING
        cutStopTimer( emptyTimer);
	emptyCounter++;
#endif


#ifdef TIMER
        ktime = cutGetTimerValue ( laptimer );
	printf(" Empty        kernel %f\n",ktime);
#endif

        //todo: correct RCount
        cutStartTimer( laptimer);

	
//        printf("Sending RCount %d\n",4*4*hBPos[i]);

#ifdef FINETUNING
        cutStartTimer( relaxTimer);
#endif

        relaxKernelPath    <<<grid, threads >>> ( 4*4*hBPos[i], dB, BCount, dBPos, dRLoc, dR, dDistR, dDist,dBucketPos, dBucketMap,lido, dSourceR, dSource);
        cudaThreadSynchronize();

#ifdef FINETUNING
        cutStopTimer( relaxTimer);
	relaxCounter++;
#endif
 

#ifdef TIMER
        ktime = cutGetTimerValue( laptimer );
        printf("Relaxing      kernel %f\n",ktime);
#endif


        cudaMemcpy( RCount,  &BCount[i],  1*sizeof(int),       cudaMemcpyDeviceToHost) ;
        cudaMemcpy( hSCount,  dSCount,  1*sizeof(int),       cudaMemcpyDeviceToHost) ;        
	

#ifdef TIMER
        printf("After RCount %d (i=%d) | SCount %d\n",RCount[0],i,hSCount[0]);
#endif


      }

//      ktime = cutGetTimerValue( laptimer );
#ifdef TIMER
      printf("Same called %d times.\n",sameCount);
#endif

//      printf("Label kernel %f\n",ktime);	
// 	printf("Before heavy label%d\n",i);      


//      cutStartTimer( laptimerrunTestGraph);

//      labelHeavyKernel    <<<grid, threads >>> ( 0, dS, hSCount[0], iw,ih, dRLoc, dR, dDistR, dDist,dBucketMap);
//      cudaThreadSynchronize();



      ktime = cutGetTimerValue( laptimer );
#ifdef TIMER
      printf("Labelling heavy kernel %f\n",ktime);
#endif

// 	printf("After heavy label%d\n",i);      
      //todo: correct SCount

#ifdef TIMER
      cutStartTimer( laptimer);
#endif
//      int temp[1];
//      cudaMemcpy( temp,  &dSCount[0],  1*sizeof(int),       cudaMemcpyDeviceToHost) ;
//      relaxKernel         <<<grid, threads >>> ( temp[0], dB, BCount, dRLoc, dR, dDistR, dDist,dBucketPos, dBucketMap,lido);

//      relaxKernel    <<<grid, threads >>> ( 4*4*hSCount[0], dB, BCount, dBPos, dRLoc, dR, dDistR, dDist,dBucketPos, dBucketMap,lido);

//      ktime = cutGetTimerValue( laptimer );
#ifdef TIMER
      printf("Relaxing heavy  kernel %f\n",ktime);
#endif

// 	printf("After heavy relax%d\n",i);      

//      printf("Done here. i = %d\n",i);


    }


    // check if kernel execution generated and error

   CUT_CHECK_ERROR("Kernel execution failed");

 

    ktime = cutGetTimerValue( nvtimer );
    cudaMemcpy( hRLoc,  dRLoc,  1024*sizeof(int),       cudaMemcpyDeviceToHost) ;    
    cudaMemcpy( hRLoc,  dRLoc,  1024*sizeof(int),       cudaMemcpyDeviceToHost) ;
    cudaMemcpy( hR,     dR,     4*1024*8*sizeof(int),   cudaMemcpyDeviceToHost) ;
    cudaMemcpy( hDistR, dDistR, 4*1024*8*sizeof(float), cudaMemcpyDeviceToHost) ;
    cudaMemcpy( hS,     dS,     1024*8*sizeof(float),   cudaMemcpyDeviceToHost) ;
    cudaMemcpy( hDist,  dDist,  GN*sizeof(int),       cudaMemcpyDeviceToHost) ;
    cudaMemcpy( hBCount, BCount,  NUMBUCKETS*sizeof(int),       cudaMemcpyDeviceToHost) ;
    cudaMemcpy( hB, dB,           2*NUMBUCKETS*sizeof(int),       cudaMemcpyDeviceToHost) ;

    cudaMemcpy( Hlido, lido, GN*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy( hR, dR, GN*sizeof(int),cudaMemcpyDeviceToHost);

    cudaMemcpy( hSource, dSource, GN*sizeof(int),  cudaMemcpyDeviceToHost);

    

  printDistances(hDist,iw,ih); 

    printPath(hSource, 10);

#ifdef FINETUNING
    printf("Label averageTime: %f ms (called %d times)\n",cutGetAverageTimerValue(labelTimer), labelCounter);
    printf("Copy  averageTime: %f ms (called %d times)\n",cutGetAverageTimerValue(copyTimer),  copyCounter);
    printf("Relax averageTime: %f ms (called %d times)\n",cutGetAverageTimerValue(relaxTimer), relaxCounter);
    printf("Empty averageTime: %f ms (called %d times)\n",cutGetAverageTimerValue(emptyTimer), emptyCounter);
    printf("Num threads: %d\n",num_threads);
    printf("Max RCount %d\n",maxRCount);
    printf("Max Distance %f\n", maxDistance);
#endif 

    for(int i=0;i<5;i++){
      printf("source[%d]= %d\n",i,hSource[i]);
    }
    for(int i=512;i<517;i++){
      printf("source[%d]= %d\n",i,hSource[i]);
    }

    printf("It took: %f ms\n", ktime);

    cudaUnbindTexture(mytex0);
    cudaUnbindTexture(mytex1);
    cudaUnbindTexture(mytex2);
    cudaUnbindTexture(mytex3);

    cudaFreeArray(array[0]);
    cudaFree(dB);
    cudaFree(dRLoc);
    cudaFree(dDist);
    cudaFree(dBucketMap);
    cudaFree(dBucketPos);

    cudaFree(BCount);

    cudaFree(dR);
    cudaFree(dDistR);
    cudaFree(dSource);
    cudaFree(dS);
    cudaFree(dSCount);
    free(dEdges[0]);
    free(dEdges[1]);
    free(dEdges[2]);
    free(dEdges[3]);
    free(nodes);
    
 
}

void empty(){
	dim3  grid( 1, 1, 1);
    dim3  threads( 256, 1, 1);
    emptyKernel    <<<grid, threads >>> ();
}

void loadEdgesFromGradient(int iw, int ih, float* gradient){

  for(int i=0;i<4;i++){
    dEdges[i] = (float*) malloc(iw*ih*sizeof(float));
    for(int j=0;j<iw*ih;j++){
      dEdges[i][j]=INF;
    }
  }

  for(int i=1;i<ih;i++){
  	for(int j=0;j<iw;j++){
  		int index = i*iw+j;
  		dEdges[RIGHT][index]    = 0.5*(1 - gradient[index]);
  		dEdges[LEFT ][index+1]  = 0.5*(1 - gradient[index+1]);
  		dEdges[DOWN] [index]    = 0.5*(1 - gradient[index]);
  		dEdges[UP]   [index-iw] = 0.5*(1 - gradient[index-iw]);
  		//printf("up %f\n",dEdges[RIGHT][index]); 
  	}
  	
  }  
    loadTexture(iw,ih,dEdges[0],array[0],&mytex0);
    loadTexture(iw,ih,dEdges[1],array[1],&mytex1);
    loadTexture(iw,ih,dEdges[2],array[2],&mytex2);
    loadTexture(iw,ih,dEdges[3],array[3],&mytex3);

}

void 
runTestGraph ( int argc, char** argv, int iw, int ih, int startNode, float* gradient){
	//loading edges from gradient
	loadEdgesFromGradient(iw,ih,gradient);	
	printf("Gradient %f\n",gradient[startNode]);
	runTest(argc,argv,iw,ih,startNode);
	printf("Delta for %d finished\n",startNode);
	
	
}

