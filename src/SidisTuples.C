#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TNtuple.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TSystemDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLeaf.h>
#include "clas12reader.h"
#include "helonline.h"
using namespace clas12;


void SetLorentzVector(TLorentzVector &p4, clas12::region_part_ptr rp, double mass){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),mass);

}

bool fillHistsDC(clas12::region_part_ptr p, TH2* dc1, TH2* dc2, TH2* dc3){
  if (p->traj(DC,DC1)->getX() != 0)
    dc1->Fill(p->traj(DC,DC1)->getX(),p->traj(DC,DC1)->getY());
  else return 0;
  if (p->traj(DC,DC2)->getX() != 0)
    dc2->Fill(p->traj(DC,DC2)->getX(),p->traj(DC,DC2)->getY());
  if (p->traj(DC,DC3)->getX() != 0)
    dc3->Fill(p->traj(DC,DC3)->getX(),p->traj(DC,DC3)->getY());
  return 1;
}


double s30 = sin(TMath::Pi()/6);
double c30 = cos(TMath::Pi()/6);
bool dcOK(clas12::region_part_ptr p){
  double x1 = p->traj(DC,DC1)->getX();
  double y1 = p->traj(DC,DC1)->getY();
  double x3 = p->traj(DC,DC3)->getX();
  double y3 = p->traj(DC,DC3)->getY();
  double cutDCa = 8; //in between sectors
  double cutDCb = 20;//inner hexagon
  if (abs(x1)<cutDCa || abs(c30*y1+s30*x1)<cutDCa || abs(c30*y1-s30*x1)<cutDCa)
    return false;
  if (abs(x1)<cutDCb && abs(c30*y1+s30*x1)<cutDCb && abs(c30*y1-s30*x1)<cutDCb)
    return false;
  if (abs(x3)<cutDCa || abs(c30*y3+s30*x3)<cutDCa || abs(c30*y3-s30*x3)<cutDCa)
    return false;
  if (abs(x3)<cutDCb && abs(c30*y3+s30*x3)<cutDCb && abs(c30*y3-s30*x3)<cutDCb)
    return false;
  //  cout << "passed" << endl;
  return true;
}

bool pcalOK(clas12::region_part_ptr p){
  double x = p->cal(PCAL)->getX();
  double y = p->cal(PCAL)->getY();

  double cutPCAL = 25; //in between sectors
  if (abs(x)<cutPCAL || abs(c30*y+s30*x)<cutPCAL || abs(c30*y-s30*x)<cutPCAL)
    return false;
  //cout << "pcal ok" <<endl;
  return true;
}

//list of cuts
double cut_ECfracmin = 0.17;
double cut_pcalmin = 0.07;
double cut_evzmin = -8;//-13;
double cut_evzmax = 1;//12;

//double cut_Wmin  = 2;
double cut_Q2min = 1;
//double cut_ymax  = 0.8;

double cut_dvzmin = -5;//-20;
double cut_dvzmax = 5;//20;

double cut_dtime_corr = 0.3;

double cut_chi2pidmin = -2.5, cut_chi2pidmax = 2.5;

//double cut_zmin = 0.3;

void toCM(TLorentzVector cm, TLorentzVector p,TLorentzVector& result){
  result = p;
  //  result.Print();
  result.RotateZ(-cm.Phi());
  result.RotateY(-cm.Theta());
  result.Boost(0,0,-cm.Beta());
  //result.Print();
  //return result;
}



void SidisTuples(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  /////////////////////////////////////
  //ignore this just getting file name!

   TString outputFile;

   TChain input("hipo");

   int maxevents=-1;
   bool useCuts = 1;
   bool isMC = 0;
   for(Int_t i=1;i<gApplication->Argc();i++){
    TString opt=gApplication->Argv(i);
    if((opt.Contains("--in="))){
      TString inputFile=opt(5,opt.Sizeof());
      cout << inputFile << endl;
      if(inputFile.EndsWith("/")){
	cout << inputFile << endl;
	inputFile=inputFile(0,inputFile.Sizeof()-1);
	
	TSystemDirectory dir(inputFile, inputFile); 
	TList *files = dir.GetListOfFiles();
	for (int j =0;j<files->GetEntries();j++){
	  TString name = files->At(j)->GetName();
	  if (name == "." || name == "..")
	    continue;
	  input.Add(inputFile+"/" +name);
	  cout << files->At(j)->GetName() <<endl;
	}
      } else
	input.Add(inputFile);
    } else if (opt.Contains("--out=")){
      outputFile = opt(6,opt.Sizeof());
    } else if (opt.Contains("--N=")){
      cout << (TString)opt(4,opt.Sizeof())<< endl;
      maxevents= ((TString)opt(4,opt.Sizeof())).Atoi();
      cout << "maxevents set to " << maxevents << endl;
    } else if (opt.Contains("--cut")){
      useCuts=1; 
    } else if (opt.EqualTo("--isMC")){
      isMC=1;
    }
   }

   clas12::clas12databases *dbc12;
   if(!isMC){
     clas12::clas12databases::SetRCDBRemoteConnection();
     dbc12 = new clas12::clas12databases();
   }
   //TString paramList = "E:helicity:e_p:e_th:e_ph:nu:Q2:x:y:W:e_DC1x:e_DC1y:e_DC2x:e_DC2y:e_DC3x:e_DC3y:e_PCALx:e_PCALy:e_vz:e_ecalfrac:e_pcal:npip:npim";
   
   //cout << paramList << endl;
   TTree* electron_tree = new TTree("electrons","electrons");//,paramList);
   //double E,ep,eth,eph,nu,Q2,x,y,W,eDC1x,eDC1y,eDC2x,eDC2y,eDC3x,eDC3y,ePCALx,ePCALy,evz;
#define leaf(name)  double name=0; tree->Branch(#name,&name,#name+(TString)"/D");
   //use this one if the variable has already been declared
#define leafx(name) tree->Branch(#name,&name,#name+(TString)"/D");
   TTree* tree = electron_tree;
   leaf(E);
   leaf(helicity);
   leaf(e_p);
   leaf(e_th);
   leaf(e_ph);
   leaf(nelectrons);
   leaf(nu);leaf(Q2);leaf(x);leaf(y);leaf(W);
   leaf(e_DC1x);leaf(e_DC2x);leaf(e_DC3x);leaf(e_DC1y);leaf(e_DC2y);leaf(e_DC3y);
   leaf(e_PCALx);leaf(e_PCALy);
   leaf(e_ecalfrac);leaf(e_pcal);
   leaf(e_vz);
   leaf(npip);leaf(npim);
   leaf(npp);leaf(npm);
   leaf(nKp);leaf(nKm);
   leaf(nh);
   leaf(z_tot);
   //cout << "made leaves" <<endl;
   
   
   TTree* hadron_tree = new TTree("hadrons","hadrons");
   tree = hadron_tree;
   leafx(nelectrons);
   leafx(E);leafx(helicity);leafx(e_p);leafx(e_th);leafx(e_ph);leafx(nu);leafx(Q2);leafx(x);leafx(y);leafx(W);
   leaf(h_chi2pid);leaf(h_pid);leaf(h_p);leaf(h_th);leaf(h_ph);leaf(h_DC1x);leaf(h_DC1y);leaf(h_DC2x);leaf(h_DC2y);leaf(h_DC3x);leaf(h_DC3y);leaf(dvz);leaf(z); leaf(h_cm_p);leaf(h_cm_th);leaf(h_cm_ph);
   leaf(h_eta); leaf(dtime); leaf(dtime_corr);
   TTree* dihadron_tree = new TTree("dihadrons","dihadrons");
   tree = dihadron_tree;
   
   leafx(E);leafx(helicity);leafx(e_p);leafx(e_th);leafx(e_ph);leafx(nu);leafx(Q2);leafx(x);leafx(y);leafx(W); 

   // macro creates fields for two hadrons
#define leaf2(name) double h1_##name=0; tree->Branch((TString)"h1_"+#name,&h1_##name,(TString)"h1_"+#name+(TString)"/D"); double h2_##name=0; tree->Branch((TString)"h2_"+#name,&h2_##name,(TString)"h2_"+#name+(TString)"/D");
   leaf2(chi2pid);leaf2(pid);leaf2(p);leaf2(th);leaf2(ph);leaf2(z);leaf2(eta);
   leaf2(cm_p);leaf2(cm_th);leaf2(cm_ph);
   leaf(pair_mass);
   leaf(diff_eta);
   leaf(diff_phi);
   leafx(nelectrons);
   leaf(diff_phi_cm);
   /*if(inputFile==TString())  {
     std::cout << " *** please provide a file name..." << std::endl;
     exit(0);
     }*/
   /////////////////////////////////////
   
   
   //cout<<"Analysing hipo file "<<inputFile<<endl;
   
   //TChain fake("hipo");
   //fake.Add(inputFile.Data());
   //get the hipo data
   //   reader.open(inputFile.Data());
   auto files=input.GetListOfFiles();
   
   //some particles
   auto db=TDatabasePDG::Instance();
   
   
   //start with dummy values
   //double recoilMass = 0;
   //by default the target mass is that of the proton.
   TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());
   TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
   TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());
   TLorentzVector cm;
   cout << "proton mass is " << db->GetParticle(2212)->Mass() << endl;
   
   
   

   gBenchmark->Start("timer");
   int count=0;
   

   for(Int_t filenum=0;filenum<files->GetEntries();filenum++){
     //create the event reader
     clas12reader c12(files->At(filenum)->GetTitle(),{0});
     if(!isMC){

       c12.connectDataBases(dbc12);
       
       clas12::rcdb_reader *rcdb = c12.rcdb();
       auto& current = rcdb->current();
       E = current.beam_energy/1000;
       TString targetName(rcdb->current().target);
       cout << "target is " << targetName << "." << endl;
       rcdb->close();
     } else{
       E=10; //Idk what the actual beam energy is for MC
     }
     
     
     TLorentzVector beam(0,0,E,E);  
     std::cout << "beam energy is " << E << std::endl;
     //  clas12reader c12(files->At(i)->GetTitle(),{0});//add tags {tag1,tag2,tag3,...}
      
      //Add some event Pid based selections
      //////////c12.AddAtLeastPid(211,1); //at least 1 pi+
      //c12.addExactPid(11,1);    //exactly 1 electron
      //c12.addExactPid(211,1);    //exactly 1 pi+
      //c12.addExactPid(-211,1);    //exactly 1 pi-
      //c12.addExactPid(2212,1);    //exactly 1 proton
      //c12.addExactPid(22,2);    //exactly 2 gamma
      //////c12.addZeroOfRestPid();  //nothing else
      //////c12.useFTBased(); //and use the Pids from RECFT

     //can also access the integrated current at this point
     //c12.scalerReader();//must call this first
     //c12.getRunBeamCharge();

     //int max = 100000;
     while(c12.next()==true){
       if(isMC){
	 auto dict = c12.getDictionary();
	 if(dict.hasSchema("MC::Event")){
	   auto schema = dict.getSchema("MC::Event");
	   hipo::bank bank(schema);
	   c12.getStructure(&bank);
	   E = bank.getFloat(schema.getEntryOrder("ebeam"),0);
	 }
	 if(dict.hasSchema("MC::Header")){
           auto schema = dict.getSchema("MC::Header");
	   hipo::bank bank(schema);
           c12.getStructure(&bank);
           helicity = bank.getFloat(schema.getEntryOrder("helicity"),0);
         }
       }
       //cout << "new event" << endl;
       //       c12.addARegionCDet();
       if(!(count %10000))
	 cout << count <<"events processed; " << maxevents << "requested"<< endl;
       count ++;

       if(count > maxevents && maxevents >0)
	 break;
	//can get an estimate of the beam current to this event
	//c12.getCurrApproxCharge();//if called c12.scalerReader();
	
        //c12.event()->getStartTime();

	
        //Loop over all particles to see how to access detector info.
	/*
 	for(auto& p : c12.getDetParticles()){
  	 //  get predefined selected information
	 p->getTime();
	 p->getDetEnergy();
	 p->getDeltaEnergy();

	 //check trigger bits
	 //	 if(c12.checkTriggerBit(25)) cout<<"MesonExTrigger"<<endl;
	 //	 else cout<<"NOT"<<endl;

	 // get any detector information (if exists for this particle)
	 // there should be a get function for any entry in the bank
	 switch(p->getRegion()) {
	 case FD :
	   p->cal(PCAL)->getEnergy();
	   p->cal(ECIN)->getEnergy();
	   p->cal(ECOUT)->getEnergy();
	   p->sci(FTOF1A)->getEnergy();
	   p->sci(FTOF1B)->getEnergy();
	   p->sci(FTOF2)->getEnergy();
	   p->trk(DC)->getSector();
	   p->che(HTCC)->getNphe();
	   p->che(LTCC)->getNphe();
	   //trajectories
	   p->traj(LTCC)->getX();
	   // p->traj(DC,DC1)->getCx();; //First layer of DC, hipo4
	   break;
	 case FT :
	   p->ft(FTCAL)->getEnergy();
	   p->ft(FTHODO)->getEnergy();
	   break;
	 case CD:
	   p->sci(CTOF)->getEnergy();
	   p->sci(CND)->getEnergy();
	   break;
	 }
	 //   covariance matrix (comment in to see!)
	 // p->covmat()->print();
	 p->cmat();
       }*/

       // get particles by type

       auto electrons=c12.getByID(11);
       //cout << "electrons" <<endl;
       //cout << electrons.size()<< "electrons" << endl;
       //if(electrons.size() > 1)
       //continue;

       //auto gammas=c12.getByID(22);
       
       //loop through all particles when searching for hadrons.
       auto parts=c12.getDetParticles();
       //cout << "parts" << endl;
       if(c12.helonline() != NULL)
	 helicity = c12.helonline()->getHelicity();
       //cout << "helicity" << endl;
       TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
       
       /*if(electrons.size() == 2){
	 cout << "electrons" << endl;
	 cout << electrons[0]->par()->getP() <<endl;
	 cout << electrons[1]->par()->getP() <<endl;
	 }*/
       nelectrons = electrons.size();
       int electrons_passCuts = 0;
       for(int i=0; i<electrons.size(); i++){
	 //if(electrons.size()>1) continue;
	 e_DC1x=electrons[i]->traj(DC,DC1)->getX();
	 e_DC1y=electrons[i]->traj(DC,DC1)->getY();
	 e_DC2x=electrons[i]->traj(DC,DC2)->getX();
	 e_DC2y=electrons[i]->traj(DC,DC2)->getY();
	 e_DC3x=electrons[i]->traj(DC,DC3)->getX();
	 e_DC3y=electrons[i]->traj(DC,DC3)->getY();
	 if(!dcOK(electrons[i]) || !pcalOK(electrons[i]))
	   continue;
	 //cout << "dc and pcal ok"<<endl;
	 //if(!dcok)
	 //  continue;
	 //cout << "electron" << endl;

	 //the electron mass is just a legend.  I could set this to zero and nothing would change.
	 SetLorentzVector(el,electrons[i], 0.000511); 
	 e_p = el.P();
	 //cout<< "e_p: " << e_p <<" "<< el.P() << endl;
	 e_th = el.Theta();
	 e_ph = el.Phi();
	 //if(theta*180/3.14159265 <=7)
	 //continue;
	 //if(p<0.01*E)
	 //continue;
	 
	 
	 double ecal = electrons[i]->getDetEnergy();
	 
	 e_ecalfrac = ecal/e_p;
	 e_pcal = electrons[i]->cal(PCAL)->getEnergy();
	 if(ecal == 0 || e_p==0 || e_pcal == 0)
	   continue;
	 e_PCALx = electrons[i]->cal(PCAL)->getX();
	 e_PCALy = electrons[i]->cal(PCAL)->getY();	

	 if(useCuts && e_ecalfrac<cut_ECfracmin)
	   continue;

	 if(e_pcal<cut_pcalmin)
	   continue;
	 
	 //cout << "/electron"<<endl;
	 
	 e_vz = electrons[i]->par()->getVz();

	 if(useCuts && (e_vz<cut_evzmin || e_vz> cut_evzmax))
	   continue;
	 
	 //done with electron id cuts


	 
	 // now for electron kinematics cuts
	 Q2 = -(beam-el)*(beam-el);
	 W = (target+beam-el).M();
	 x = Q2/(2*target.M()*(beam.E()-el.E()));
	 nu = (beam.E()-el.E());
	 y = nu/E;
	 
	 if(useCuts && Q2<cut_Q2min)
           continue;
	 
	 //if(useCuts && y>cut_ymax)
	 //continue;
	 cm = beam+target-el;
	 //cout << "cm:"<<endl;
	 //cm.Print();
	 npip = 0;
	 npim = 0;
	 npp = 0;
	 npm = 0;
	 nKp = 0;
	 nKm = 0;
	 nh = 0;
	 TLorentzVector had;
	 TLorentzVector h1;
	 TLorentzVector h2;

	 z_tot = 0;
	 bool found_leader = 0, found_second=0;
	 //loop through all particles, only choosing charged hadrons. 
	 for(int j =0; j<parts.size();j++){

	   auto h = parts[j];
	   h_pid = h->getPid();
	   if(h_pid != 211 && h_pid != -211 && h_pid != 2212 && h_pid != -2212 && h_pid != 321 && h_pid != -321)
	     continue;
	   h_DC1x=h->traj(DC,DC1)->getX();
	   h_DC1y=h->traj(DC,DC1)->getY();
	   h_DC2x=h->traj(DC,DC2)->getX();
	   h_DC2y=h->traj(DC,DC2)->getY();
	   h_DC3x=h->traj(DC,DC3)->getX();
	   h_DC3y=h->traj(DC,DC3)->getY();
	   
	   dtime = electrons[i]->getTime()-h->getTime();
	   
	   //cout << "dc"<<endl;
	   if(!dcOK(h))
	     continue;

           //bool dcok = fillHistsDC(pips[j],hpipdc1xy,hpipdc2xy,hpipdc3xy);
	   //if(!dcok)
	   //continue;
	   //cout << "pip" << endl;
	   double mass = db->GetParticle(h_pid)->Mass();
	   SetLorentzVector(had,h, mass);
	   double c = 29.9792458; //cm/ns 
	   dtime_corr =dtime-electrons[i]->getPath()/c+h->getPath()/(had.Beta()*c);
	   if(abs(dtime_corr) > cut_dtime_corr)
	     continue;
	   h_p = had.P();
	   h_th = had.Theta();
	   h_eta = had.PseudoRapidity();
	   h_ph = had.Phi();
	   
	   //h_pid = h->par()->getPid();
	   h_chi2pid = h->par()->getChi2Pid();
	   //cout << "chi2pid" << endl;
	   if(useCuts && abs(h_chi2pid) > cut_chi2pidmax)
	     continue;
	   
	   dvz = electrons[i]->par()->getVz()-h->par()->getVz();
	   
	   if(useCuts && (dvz < cut_dvzmin || dvz > cut_dvzmax))
	     continue;
	   
	   z = had.E()/nu;
	   
	   TLorentzVector h_cm;
	   toCM(cm, had,h_cm);
	   h_cm_p = h_cm.P();
	   h_cm_th = h_cm.Theta();
	   h_cm_ph = h_cm.Phi();
	   //cout << h_cm_p << " " << pi_cm_th << " " << pi_cm_ph << endl;
	   hadron_tree->Fill();
	   
	   //leading pion in a high-z pion, and a second hadron of any type
	   if(z > 0.5 && abs(h_pid)==211){
	     h1_pid = h_pid;
	     h1_p = h_p;
	     h1_th = h_th;
	     h1_eta = h_eta;
	     h1_ph = h_ph;
	     h1_chi2pid = h_chi2pid;
	     
	     h1_cm_p = h_cm_p;
	     h1_cm_th = h_cm_th;
	     h1_cm_ph =h_cm_ph;
	     h1_z = z;
	     h1 = had;
	     found_leader = 1;
	   } 
	   else {
	     h2_pid = h_pid;
	     h2_p = h_p;
	     h2_th = h_th;
	     h2_eta = h_eta;
	     h2_ph = h_ph;
	     h2_chi2pid = h_chi2pid;

	     h2_cm_p = h_cm_p;
	     h2_cm_th =h_cm_th;
	     h2_cm_ph =h_cm_ph;
	     h2_z = z;
	     h2 =had;
	     found_second = 1;
	   }
	   if(found_leader && found_second){
	     //cout << "masses "  << p1.M() << "  " << p2.M() <<endl;
	     pair_mass = (h1+h2).M();
	     diff_eta = h1_eta-h2_eta;
	     double PI = TMath::Pi();
	     diff_phi = h1_ph-h2_ph;
	     if(diff_phi<-PI)
	       diff_phi+=2*PI;
	     if(diff_phi>PI)
	       diff_phi-=2*PI;
	     
	     diff_phi_cm = h1_cm_ph-h2_cm_ph;
	     if(diff_phi_cm<-PI)
	       diff_phi_cm+=2*PI;
	     if(diff_phi_cm>PI)
	       diff_phi_cm-=2*PI;
	     
	     dihadron_tree->Fill();
	   } 
	   
	   if(h_pid == 211)
	     npip++;
	   else if(h_pid == -211)
	     npim++;
	   else if(h_pid == 2212)
             npp++;
	   else if(h_pid == -2212)
             npm++;
	   else if(h_pid == 321)
             nKp++;
           else if(h_pid == -321)
             nKm++;
	   nh++;
	   z_tot+=z;
	 }
	 //cout<< "e_p: " <<e_p<<endl;
	 electron_tree->Fill();
	 electrons_passCuts++;
       }
       //cout << electrons_passCuts << " electrons pass cuts" << endl;
       //cout << "end of file"<< endl;

       
     } 
      
   }
   gBenchmark->Stop("timer");
   gBenchmark->Print("timer");
   //hm2gCut->SetLineColor(2);
   //hm2gCut->DrawCopy("same");
   TFile *f = new TFile(outputFile,"RECREATE");
   electron_tree->Write();
   hadron_tree->Write();
   dihadron_tree->Write();
   f->Close();
   
   
   

   auto finish = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = finish - start;
   std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<count<< " s\n";

}
