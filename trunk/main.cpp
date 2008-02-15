#include <QApplication>
#include <QPushButton>
#include "Window.h"

void
runTest( int argc, char** argv,int iw, int ih, int startNode);

void 
runTestGraph( int argc, char** argv,int iw, int ih, int startNode, float* gradient);

void empty();

int main(int argc, char* argv[]){
	//only to initialize the environment
	empty();
	
	/*float* dummy;
  if(argc <2){
    printf("Usage: program_name graph_file graph_width graph_height start_node\n");
    return 0;
  }*/
  //printf( default iw and ih is 512
/*  if(argc == 2){
    runTest( argc, argv,512,512,0);
  }
  else{
    int iw, ih;
    sscanf(argv[2],"%d",&iw);
    sscanf(argv[3],"%d",&ih);      
    int startNode = 0;
    if(argc >= 5)
      sscanf(argv[4],"%d",&startNode);
    //runTestGraph( argc, argv,iw,ih,startNode,dummy);
  }*/  

  QApplication app(argc,argv);      
  Window window;
  window.show();
  
  return app.exec();
}
