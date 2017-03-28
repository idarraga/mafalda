#include <TControlBar.h>
#include <TApplication.h>
#include <TROOT.h>
#include <iostream>

class BControl : public TObject {
   
protected:
   TControlBar   *fControlBar;  // control bar
   
public:
   BControl();
   //virtual ~BControl();
   ~BControl();

   void Build();
   void PrintEventStats();
   void seekForward();
   void seekBack();
   void ShowResults();
   void jump100F();
   void jump100B();

   TApplication * m_theApp;

   ClassDef(BControl,0)
};

