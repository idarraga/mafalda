#include "control.h"

extern TApplication * g_theApp;
extern Int_t g_direction;

ClassImp(BControl)

BControl::BControl()
{
   fControlBar = 0;
   ShowResults();
}

BControl::~BControl(){
   if (fControlBar) delete fControlBar;
}

void BControl::Build(){

   fControlBar = new TControlBar("vertical", "MPX Viewer"); // Orientation & title.
   
   fControlBar->AddButton("seek >>", "BControl::seekForward()", "Forward");
   fControlBar->AddButton("<<", "BControl::seekBack()", "Back");
   fControlBar->AddButton(">> 100", "BControl::jump100F()", "Jump 100");
   fControlBar->AddButton("<< 100", "BControl::jump100B()", "Jump 100");
   fControlBar->AddButton("PrintEventStats", "BControl::PrintEventStats()", "Calls PrintEventStats");
   fControlBar->AddButton("Quit", ".q", "Quits ROOT");
   fControlBar->Show();
}

void BControl::PrintEventStats(){

   printf("In PrintEventStats()\n");

}

void BControl::jump100F(){
  g_theApp->Terminate();
  g_direction = 51;
}
void BControl::jump100B(){
  g_theApp->Terminate();
  g_direction = 50;
}
void BControl::seekForward(){

  g_theApp->Terminate();
  g_direction = 1;

}

void BControl::seekBack(){

  g_theApp->Terminate();
  g_direction = 0;

}

void BControl::ShowResults(){

   Build();
   printf("In ShowResults() buid the Control Bar\n");
}

void cbar1() {
   new BControl();
}
