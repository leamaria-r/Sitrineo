void sitrineocumul(Int_t runNumber=1210, Int_t nEvt=100000) {
  // Macro to start continuous SITRINEO analysis with evt by evt display
  // launch with: nohup TAF code/macros/sitrineocumul.C &
  //
  // JB 2021/10/02

  gTAF->InitSession(runNumber,1, 0, "./config/sitrineo-m2-woB-cp.cfg");
//  gTAF->SetDebug(-2);

  gTAF->GetRaw()->SitrineoCumul(nEvt, 2);

  gApplication->Terminate();
}
