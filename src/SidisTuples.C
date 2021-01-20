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
//#include "Math/Vector4Dfwd.h"
using namespace clas12;

int debug = 0;
double defval = 0;

void SetLorentzVector(TLorentzVector &p4, clas12::region_part_ptr rp, double mass){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),mass);

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
double cut_Wmin = 2;
double cut_ymax  = 0.85;

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

double angle(double phi){
  double pi = TMath::Pi();
  while(phi>pi)
    phi-=2*pi;
  while(phi<-pi)
    phi+=2*pi;
  return phi;
}


void SidisTuples(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  /////////////////////////////////////
  //ignore this just getting file name!

   TString outputFile;
   TChain input("hipo");
   int skipEvents=0;
   int maxevents=-1;
   bool useCuts = 1;
   bool isMC = 0;
   bool createDipionTree = 0;
   bool createElectronTree = 0;
   bool createDihadronTree = 0;
   bool createHadronTree = 0;
   for(Int_t ii=1;ii<gApplication->Argc();ii++){
    TString opt=gApplication->Argv(ii);
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
    } else if(opt.Contains("--skipEvents=")){
      skipEvents = ((TString)opt(13,opt.Sizeof())).Atoi();
      cout << "skipping the first " << skipEvents << " events" << endl;
    } else if (opt.Contains("--cut")){
      useCuts=1; 
    } else if (opt.EqualTo("--isMC")){
      isMC=1;
    } else if (opt.Contains("--includeDipions")){
      createDipionTree = 1;
      cout << "create dipion tree" <<endl;
    } else if (opt.Contains("--includeDihadrons")){
	createDihadronTree = 1;
	cout << "create dihadron tree" <<endl;
    } else if (opt.Contains("--includeHadrons")){
      createHadronTree = 1;
      cout << "create hadron tree" <<endl;
    } else if (opt.Contains("--includeElectrons")){
      createElectronTree = 1;
      cout << "create electron tree" <<endl;
    } else if (opt.Contains("--debug")){
      debug=1;
    }

   }

   if(!(createDihadronTree||createHadronTree||createDipionTree||createElectronTree)){
     cout << "no trees selected!  Aborting" << endl;
     exit(0);
   }

   clas12::clas12databases *dbc12;
   if(!isMC){
     clas12::clas12databases::SetRCDBRemoteConnection();
     dbc12 = new clas12::clas12databases();
   }
   //macro for declaring a variable and adding it to the current tree
#define leaf(name)  double name=0; if(tree != NULL) tree->Branch(#name,&name,#name+(TString)"/D");
   //use this one if the variable has already been declared
#define leafx(name) if(tree != NULL) tree->Branch(#name,&name,#name+(TString)"/D");
   TTree* tree = NULL;

   TTree* electron_tree= createElectronTree ? new TTree("electrons","electrons") : NULL;//,paramList);
   tree = electron_tree;
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
   leaf(ntracks);
   leaf(z_tot);

   
   // don't add these variables to the tree unless running MC 
   if(!isMC)
     tree = NULL;
   leaf(e_truth_pid);leaf(e_truth_p);leaf(e_truth_th);leaf(e_truth_ph);leaf(e_truth_vx);leaf(e_truth_vy); leaf(e_truth_vz);
   

   //cout << "made leaves" <<endl;
   
   
   TTree* hadron_tree = createHadronTree ? new TTree("hadrons","hadrons") : NULL;
   tree = hadron_tree;
   leafx(nelectrons);
   leafx(E);leafx(helicity);leafx(e_p);leafx(e_th);leafx(e_ph);leafx(nu);leafx(Q2);leafx(x);leafx(y);leafx(W);leafx(ntracks);
   leaf(h_chi2pid);leaf(h_pid);leaf(h_p);leaf(h_th);leaf(h_ph);leaf(h_DC1x);leaf(h_DC1y);leaf(h_DC2x);leaf(h_DC2y);leaf(h_DC3x);leaf(h_DC3y);leaf(dvz);leaf(z); leaf(h_cm_p);leaf(h_cm_th);leaf(h_cm_ph);leaf(h_cm_eta);leaf(h_cm_pt);
   leaf(h_eta); leaf(dtime); leaf(dtime_corr); leaf(missing_mass);

   
   if(!isMC) tree = NULL;
   leafx(e_truth_pid);leafx(e_truth_p);leafx(e_truth_th);leafx(e_truth_ph);
   leaf(h_truth_pid);leaf(h_truth_p);leaf(h_truth_th);leaf(h_truth_ph);leaf(h_truth_cm_p);leaf(h_truth_cm_th);leaf(h_truth_cm_ph);leaf(h_truth_cm_eta);leaf(h_truth_cm_pt);leaf(h_truth_z);leaf(missing_mass_truth);
   

   TTree* dihadron_tree = createDihadronTree ? new TTree("dihadrons","dihadrons") : NULL;
   tree = dihadron_tree;   
   leafx(E);leafx(helicity);leafx(e_p);leafx(e_th);leafx(e_ph);leafx(nu);leafx(Q2);leafx(x);leafx(y);leafx(W); leafx(ntracks); 
   

   // macro creates fields for two hadrons
#define leaf2(name) double h1_##name=0; double h2_##name=0; if(tree != NULL) {tree->Branch((TString)"h1_"+#name,&h1_##name,(TString)"h1_"+#name+(TString)"/D"); tree->Branch((TString)"h2_"+#name,&h2_##name,(TString)"h2_"+#name+(TString)"/D");}
   leaf2(chi2pid);leaf2(pid);leaf2(p);leaf2(th);leaf2(ph);leaf2(z);leaf2(eta);
   leaf2(cm_p);leaf2(cm_th);leaf2(cm_ph);leaf2(cm_eta);leaf2(cm_pt);
   leaf(pair_mass);leaf(mx_eh1h2x);leaf(mx_eh1x);leaf(mx_eh2x);
   leaf(diff_eta);
   leaf(diff_phi);
   leafx(nelectrons);
   leaf(diff_phi_cm);
   leaf(diff_eta_cm);



   
   if(!isMC) tree = NULL;
   leafx(e_truth_pid);leafx(e_truth_p);leafx(e_truth_th);leafx(e_truth_ph);
   leaf(h1_truth_z);
   leaf(h1_truth_pid);leaf(h1_truth_p);leaf(h1_truth_th);leaf(h1_truth_ph);leaf(h1_truth_cm_p);
   leaf(h1_truth_cm_th);leaf(h1_truth_cm_ph);leaf(h1_truth_cm_eta);leaf(h1_truth_cm_pt);
   leaf(h2_truth_pid);leaf(h2_truth_p);leaf(h2_truth_th);leaf(h2_truth_ph);leaf(h2_truth_cm_p);leaf(h2_truth_z);
   leaf(h2_truth_cm_th);leaf(h2_truth_cm_ph);leaf(h2_truth_cm_eta);leaf(h2_truth_cm_pt);
   leaf(diff_phi_cm_truth);leaf(diff_eta_cm_truth);leaf(pair_mass_truth);leaf(mx_eh1h2x_truth);leaf(mx_eh1x_truth);leaf(mx_eh2x_truth);
   

   // A small tree for storing dipion events, without requiring one of them to be leading
   //
     
   TTree* dipion_tree = createDipionTree ? new TTree("dipions","dipions") : NULL;
   tree = dipion_tree;
   leafx(nelectrons);
   leafx(E);leafx(helicity);leafx(e_p);leafx(e_th);leafx(e_ph);leafx(nu);leafx(Q2);leafx(x);leafx(y);leafx(W);leafx(ntracks);
   leaf(pi1_p);leaf(pi2_p);leaf(pi1_th);leaf(pi2_th);leaf(pi1_ph);leaf(pi2_ph);
   leaf(pi1_cm_eta);leaf(pi2_cm_eta);leaf(pi1_cm_pt);leaf(pi2_cm_pt);leaf(pi1_z);leaf(pi2_z);leaf(pi1_pid);leaf(pi2_pid);leaf(pi1_cm_ph);leaf(pi2_cm_ph);

   leafx(diff_phi_cm);leafx(diff_eta_cm);leafx(pair_mass);leaf(mx_epi1pi2x);leaf(mx_epi1x);leaf(mx_epi2x);
   
   if(!isMC) tree = NULL;
   leafx(e_truth_pid);leafx(e_truth_p);leafx(e_truth_th);leafx(e_truth_ph);
   
   leaf(pi1_truth_pid);leaf(pi1_truth_p);leaf(pi1_truth_th);leaf(pi1_truth_ph);leaf(pi1_truth_z);
   leaf(pi1_truth_cm_ph);leaf(pi1_truth_cm_eta);leaf(pi1_truth_cm_pt);
   leaf(pi2_truth_pid);leaf(pi2_truth_p);leaf(pi2_truth_th);leaf(pi2_truth_ph);leaf(pi2_truth_z);
   leaf(pi2_truth_cm_ph);leaf(pi2_truth_cm_eta);leaf(pi2_truth_cm_pt);
   leafx(diff_phi_cm_truth);leafx(diff_eta_cm_truth);leafx(pair_mass_truth);leaf(mx_epi1pi2x_truth);leaf(mx_epi1x_truth);leaf(mx_epi2x_truth);
   
   

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
       std::cout << "beam energy is " << E << std::endl;
       rcdb->close();
     } else{
       E=0; //get the beam energy later
     }
     
     
     TLorentzVector beam(0,0,E,E);  

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


     // for data, this step can be done.
     // with MC, there are events with 
     // electrons that fail recon,
     // so these must not be excluded
     if(!isMC)
       c12.addAtLeastPid(11,1);
     while(c12.next()==true){
       if(isMC && E==0){
	 auto dict = c12.getDictionary();
	 if(dict.hasSchema("MC::Event")){
	   auto schema = dict.getSchema("MC::Event");
	   hipo::bank bank(schema);
	   c12.getStructure(&bank);
	   E = bank.getFloat(schema.getEntryOrder("ebeam"),0);
	   beam = {0, 0, E, E};
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
       if(!((count-skipEvents) %10000) && count >= skipEvents)
	 cout << (count-skipEvents) <<"events processed; " << maxevents << "requested"<< endl;
       //count ++;

       if(count-skipEvents > maxevents && maxevents >0)
	 break;
       if(count < skipEvents){
	 count ++;
	 continue;
       }
       count ++;
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
       
       if(c12.helonline() != NULL)
	 helicity = c12.helonline()->getHelicity();
       //if(debug) cout << "helicity" << endl;
       TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
       
       mcpar_ptr mcparts;
       if(isMC){
	 mcparts = c12.mcparts();
       }
       auto parts=c12.getDetParticles(); 
       /*if(electrons.size() == 2){
	 if(debug) cout << "electrons" << endl;
	 if(debug) cout << electrons[0]->par()->getP() <<endl;
	 if(debug) cout << electrons[1]->par()->getP() <<endl;
	 }*/
       nelectrons = electrons.size();
       int electrons_passCuts = 0;
       vector<int> matchedMCindices = {};
       //if(debug) cout <<"CHECK 0: " << nelectrons << " electrons" << endl;
       for(int i=0; i<nelectrons; i++){
	 if(debug) cout << "starting electron loop" << endl;
	 //if(debug) cout << "CHECK 0.5" << endl;
	 //if(electrons.size()>1) continue;
	 e_DC1x=electrons[i]->traj(DC,DC1)->getX();
	 e_DC1y=electrons[i]->traj(DC,DC1)->getY();
	 e_DC2x=electrons[i]->traj(DC,DC2)->getX();
	 e_DC2y=electrons[i]->traj(DC,DC2)->getY();
	 e_DC3x=electrons[i]->traj(DC,DC3)->getX();
	 e_DC3y=electrons[i]->traj(DC,DC3)->getY();
	 if(!dcOK(electrons[i]) || !pcalOK(electrons[i]))
	   continue;
	 //if(debug) cout << "dc and pcal ok"<<endl;
	 //if(!dcok)
	 //  continue;
	 //if(debug) cout << "electron" << endl;

	 //the electron mass is a myth.  I could set this to zero and nothing would change.
	 //if(debug) cout << "CHECK 1" <<endl;
	 SetLorentzVector(el,electrons[i], 0.000511); 
	 e_p = el.P();
	 //if(debug) cout<< "e_p: " << e_p <<" "<< el.P() << endl;
	 e_th = el.Theta();
	 e_ph = el.Phi();
	 //if(theta*180/3.14159265 <=7)
	 //continue;
	 //if(p<0.01*E)
	 //continue;
	 
	 //if(debug) cout << "CHECK 2"<<endl;
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
	 
	 //if(debug) cout << "/electron"<<endl;
	 
	 e_vz = electrons[i]->par()->getVz();

	 if(useCuts && (e_vz<cut_evzmin || e_vz> cut_evzmax))
	   continue;
	 
	 //done with electron id cuts
	 if(debug) cout << "check 2" <<endl;

	 //if(debug) cout << "now for electron kinematics cuts"<<endl;
	 // now for electron kinematics cuts
	 Q2 = -(beam-el)*(beam-el);
	 W = (target+beam-el).M();
	 x = Q2/(2*target.M()*(beam.E()-el.E()));
	 nu = (beam.E()-el.E());
	 y = nu/E;
	 
	 if(useCuts && (Q2<cut_Q2min || W<cut_Wmin || y > cut_ymax))
           continue;
	 
	 //if(useCuts && y>cut_ymax)
	 //continue;
	 cm = beam+target-el;

	 
	 auto parts=c12.getDetParticles();
	 ntracks=0;
	 if(debug) cout << "counting tracks"<<endl;
	 for(int kk = 0; kk<parts.size();kk++){
	   auto part = parts[kk];
	   int pid = part->getPid();

	   if(db->GetParticle(pid) == NULL || db->GetParticle(pid)->Charge() == 0)
	     continue;
	   double dtime = electrons[i]->getTime()-part->getTime();
	   
	   double mass = db->GetParticle(pid)->Mass();
	   double c = 29.9792458; //cm/ns
	   double dtime_corr =dtime-electrons[i]->getPath()/c+part->getPath()/(part->getBetaFromP()*c);
	   if(abs(dtime_corr) > cut_dtime_corr)
	     continue;
	   double dvz = electrons[i]->par()->getVz()-part->par()->getVz();
	   if(dvz < cut_dvzmin || dvz > cut_dvzmax)
	     continue;
	   ntracks++;
	   
	   
	 }
	 if(debug) cout << ntracks << " tracks" << endl;
	 //if(debug) cout << "parts" << endl;
	 


	 //if(debug) cout << "cm:"<<endl;
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
	 TLorentzVector e_truth;
	 TLorentzVector cm_truth;

	 if(isMC){
	   if(debug) cout << "filling electron MC" << endl;
	   double best_match_diff = 99999;
	   e_truth_p = 0;
	   e_truth_th = 0;
	   e_truth_ph = 0;
	   int kbest = -1;
	   TVector3 mc;
	   for(int k = 0; k<mcparts->getRows();k++){
	     //if(mcparts->getPid(k) != 11)
	     //continue;
	     //if(debug) cout << "e passed charge"<<endl;
	     mc = {mcparts->getPx(k),mcparts->getPy(k),mcparts->getPz(k)};
	     //if(debug) cout << mc.Theta() << "\t" << e_th << endl;
	     if(abs(mc.Theta()-e_th)>1*TMath::Pi()/180){
	       continue;
	     }
	     //if(debug) cout << "e passed phi" <<endl;
	     if(abs(mc.Theta()-e_th)>1*TMath::Pi()/180)
	       continue;
	     //if(debug) cout << "e passed phi" <<endl;
	     double diff = hypot(angle(mc.Phi()-e_ph)*sin(e_th),mc.Theta()-e_th);
	     //closest match which has negative charge
	     if(diff < best_match_diff){
	       kbest = k;
	       best_match_diff = diff;
	     }
	   }
	   if(kbest >= 0){
	     matchedMCindices.push_back(kbest);
	     // the electron mass is a myth.  
	     // I could set it to zero and nothing would change in the analysis
	     e_truth.SetXYZM(mcparts->getPx(kbest),mcparts->getPy(kbest),mcparts->getPz(kbest),0.000511);
	     e_truth_pid = mcparts->getPid(kbest);
	     e_truth_p = e_truth.P();
	     e_truth_th = e_truth.Theta();
	     e_truth_ph = e_truth.Phi();
	     cm_truth = beam+target-e_truth;
	     e_truth_vx = mcparts->getVx(kbest);
	     e_truth_vy = mcparts->getVy(kbest);
	     e_truth_vz = mcparts->getVz(kbest);
	   }
	 }
	 if(debug) cout << " check 1"<< endl;

	 
	 bool found_leader = 0, found_second=0;
	 //loop through all particles, only choosing charged hadrons.
	 if(hadron_tree != NULL || dihadron_tree != NULL || dipion_tree != NULL){
	 for(int j =0; j<parts.size();j++){
	   if(debug) cout << "hadrons loop"<< endl;
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
	   
	   //if(debug) cout << "dc"<<endl;
	   if(!dcOK(h))
	     continue;

           //bool dcok = fillHistsDC(pips[j],hpipdc1xy,hpipdc2xy,hpipdc3xy);
	   //if(!dcok)
	   //continue;
	   //if(debug) cout << "pip" << endl;
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
	   //if(debug) cout << "chi2pid" << endl;
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
	   h_cm_eta = h_cm.PseudoRapidity();
	   h_cm_ph = h_cm.Phi();
	   h_cm_pt = h_cm.Pt();
	   //if(debug) cout << h_cm_p << " " << pi_cm_th << " " << pi_cm_ph << endl;
	   
	   TLorentzVector h_truth;
	   if(isMC){
	     double best_match_diff = 9999;
	     h_truth_pid = 0;
	     h_truth_p = 0;
	     h_truth_th = 0;
	     h_truth_ph = 0;
	     h_truth_cm_p = 0;
	     h_truth_cm_th = 0;
	     h_truth_cm_eta = 0;
	     h_truth_cm_ph = 0;
	     h_truth_cm_pt = 0;
	     int kbest = -1;
	     for(int k = 0; k<mcparts->getRows();k++){
	       //if(debug) cout << "hadron" <<endl;
	       if(db->GetParticle(h_pid)->Charge() != db->GetParticle(mcparts->getPid(k))->Charge())
	       continue;
	       //if(debug) cout << "passed charge" << endl;
	       TVector3 mc(mcparts->getPx(k),mcparts->getPy(k),mcparts->getPz(k));
	       if(abs(mc.Theta()-h_th)>1*TMath::Pi()/180){
		 continue;
	       }
	       //if(debug) cout << "passed theta" << endl;
	       if(abs(angle(mc.Phi()-h_ph))>3*TMath::Pi()/180){
		 continue;
	       }
	       //if(debug) cout << "passed phi" <<endl;
	       double diff = hypot(angle(mc.Phi()-h_ph)*sin(h_th),mc.Theta()-h_th);
	       //closest match which has negative charge
	       //if(debug) cout << "diff " << diff << endl;
	       if(diff < best_match_diff){
		 kbest = k;
		 best_match_diff = diff;

	       }
	     }
	     if(kbest >= 0){
	       matchedMCindices.push_back(kbest);
	       TLorentzVector h_truth_cm;
	       
	       h_truth.SetXYZM(mcparts->getPx(kbest),mcparts->getPy(kbest),mcparts->getPz(kbest),db->GetParticle(mcparts->getPid(kbest))->Mass());
	       toCM(cm_truth, h_truth,h_truth_cm);
	       h_truth_pid = mcparts->getPid(kbest);
	       h_truth_p = h_truth.P();
	       h_truth_th = h_truth.Theta();
	       h_truth_ph = h_truth.Phi();
	       h_truth_z = h_truth_p/(E-e_truth_p);

	       h_truth_cm_p = h_truth_cm.P();
	       h_truth_cm_th = h_truth_cm.Theta();
	       h_truth_cm_eta = h_truth_cm.PseudoRapidity();
	       h_truth_cm_ph = h_truth_cm.Phi();
	       h_truth_cm_pt = h_truth_cm.Pt();
	     } 
	   }
	   if(hadron_tree != NULL){
	     missing_mass = (beam+target-el-had).M();
	     hadron_tree->Fill();
	   }
	   
	   
	   if(debug) cout << "filled hadron tree"<< endl;
	   //leading pion is a high-z pion, and a second hadron of any type
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
	     h1_cm_eta = h_cm_eta;
	     h1_cm_pt = h_cm_pt;
	     h1_z = z;
	     h1 = had;
	     found_leader = 1;
	     if(isMC){
	       h1_truth_pid = h_truth_pid;
               h1_truth_p = h_truth_p;
               h1_truth_th = h_truth_th;
               h1_truth_ph = h_truth_ph;
	       h1_truth_z = h_truth_z;

	       h1_truth_cm_p = h_truth_cm_p;
               h1_truth_cm_th = h_truth_cm_th;
               h1_truth_cm_eta = h_truth_cm_eta;
               h1_truth_cm_ph = h_truth_cm_ph;
               h1_truth_cm_pt = h_truth_cm_pt;
	     }
	   } 
	   else {
	     h2_pid = h_pid;
	     h2_p = h_p;
	     h2_th = h_th;
	     h2_eta = h_eta;
	     h2_ph = h_ph;
	     h2_chi2pid = h_chi2pid;
	     h2_truth_z = h_truth_z;

	     h2_cm_p = h_cm_p;
	     h2_cm_th =h_cm_th;
	     h2_cm_ph =h_cm_ph;
	     h2_cm_eta = h_cm_eta;
	     h2_cm_pt = h_cm_pt;
	     h2_z = z;
	     h2 =had;
	     found_second = 1;
	     if(isMC){
               h2_truth_pid = h_truth_pid;
               h2_truth_p = h_truth_p;
               h2_truth_th = h_truth_th;
               h2_truth_ph = h_truth_ph;

               h2_truth_cm_p = h_truth_cm_p;
               h2_truth_cm_th = h_truth_cm_th;
               h2_truth_cm_eta = h_truth_cm_eta;
               h2_truth_cm_ph = h_truth_cm_ph;
               h2_truth_cm_pt = h_truth_cm_pt;
             }
	   }
	   if(found_leader && found_second){
	     //if(debug) cout << "masses "  << p1.M() << "  " << p2.M() <<endl;
	     pair_mass = (h1+h2).M();
	     mx_eh1h2x = (beam+target-el-h1-h2).M();
	     mx_eh1x = (beam+target-el-h1).M();
	     mx_eh2x = (beam+target-el-h2).M();
	     diff_eta = h1_eta-h2_eta;
	     double PI = TMath::Pi();
	     diff_phi = h1_ph-h2_ph;
	     if(diff_phi<-PI)
	       diff_phi+=2*PI;
	     if(diff_phi>PI)
	       diff_phi-=2*PI;
	     
	     diff_phi_cm = h2_cm_ph-h1_cm_ph;
	     if(diff_phi_cm<-PI)
	       diff_phi_cm+=2*PI;
	     if(diff_phi_cm>PI)
	       diff_phi_cm-=2*PI;
	     diff_eta_cm = h2_cm_eta-h1_cm_eta;
	     if(isMC){
	       diff_phi_cm_truth = h2_truth_cm_ph-h1_truth_cm_ph;
	       diff_eta_cm_truth = h2_truth_cm_eta-h1_truth_cm_eta;
	     }
	     
	     if(dihadron_tree != NULL)
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
	 if(electron_tree != NULL)
	   electron_tree->Fill();
	 electrons_passCuts++;

	 //dipion tree
	 if(debug) cout << "checking for dipions" << endl;
	 if(abs(h_pid) == 211 && createDipionTree){
	   for(int k = 0; k<parts.size();k++){
	     if(debug) cout << "check1" << endl;
	     auto h2 = parts[k];
	     int pid = h2->getPid();
	     TLorentzVector had2;
	     TLorentzVector had2_cm;
	     SetLorentzVector(had2,h2, db->GetParticle(211)->Mass());
	     if(abs(pid) == 211 && had2.E()/nu < z){
	       pi1_z = z;
	       pi1_p = had.P();
	       pi1_th = had.Theta();
	       pi1_ph = had.Phi();
	       toCM(cm, had2,had2_cm);
	       pi1_cm_pt = h_cm_pt;
	       pi1_cm_ph = h_cm_ph;
	       pi1_cm_eta = h_cm_eta;
	       pi1_pid = h_pid;
	       
	       pi2_z = had2.E()/nu;
	       pi2_p = had2.P();
	       pi2_th = had2.Theta();
	       pi2_ph = had2.Phi();
	       pi2_cm_pt = had2_cm.Pt();
	       pi2_cm_ph = had2_cm.Phi();
	       pi2_cm_eta = had2_cm.Eta();
	       pi2_pid = pid;
	       diff_phi_cm = angle(pi1_cm_ph-pi2_cm_ph);
	       diff_eta_cm = pi1_cm_eta-pi2_cm_eta;
	       pair_mass = (had2+had).M();
	       mx_epi1pi2x = (beam+target-el-had-had2).M();
	       mx_epi1x = (beam+target-el-had).M();
	       mx_epi2x = (beam+target-el-had2).M();
	       if(isMC){
		 pi1_truth_z=0;
		 pi1_truth_pid=0;
		 pi1_truth_p=0;
		 pi1_truth_th=0;
		 pi1_truth_ph=0;
		 pi1_truth_cm_pt=0;
		 pi1_truth_cm_ph=0;
		 pi1_truth_cm_eta=0;

		 pi2_truth_z=0;
                 pi2_truth_pid=0;
                 pi2_truth_p=0;
                 pi2_truth_th=0;
                 pi2_truth_ph=0;
                 pi2_truth_cm_pt=0;
                 pi2_truth_cm_ph=0;
                 pi2_truth_cm_eta=0;

		 double p1_p = had.P();
		 double pi2_th = had2.Theta();
		 double pi2_ph = had2.Phi();
		 double best_match_diff = 9999;
		 int kbest1 = -1;
		 int kbest2 = -1;
		 for(int k = 0; k<mcparts->getRows();k++){
		     //if(debug) cout << "hadron" <<endl;
		   TVector3 mc(mcparts->getPx(k),mcparts->getPy(k),mcparts->getPz(k));
		   if(abs(mc.Theta()-pi2_th)>1*TMath::Pi()/180){
		     continue;
		   }
		   if(abs(angle(mc.Phi()-pi2_ph))>3*TMath::Pi()/180){
		     continue;
		   }
		   //if(debug) cout << "passed phi" <<endl;                                   
		   double diff = hypot(angle(mc.Phi()-pi2_ph)*sin(pi2_th),mc.Theta()-pi2_th);
		   //closest match which has negative charge                        
		   //if(debug) cout << "diff " << diff << endl;                               
		   if(diff < best_match_diff){
		     kbest1 = k;
		     best_match_diff = diff;
		   }
		 }
		 if(debug) cout << "loop on mc" << endl;
		 best_match_diff = 9999;
		 for(int k = 0; k<mcparts->getRows();k++){
		   if(k == kbest1)
		     continue;
		   //if(debug) cout << "hadron" <<endl;                                                                  
                   TVector3 mc(mcparts->getPx(k),mcparts->getPy(k),mcparts->getPz(k));
                   if(abs(mc.Theta()-h_th)>1*TMath::Pi()/180){
                     continue;
                   }
                   if(abs(angle(mc.Phi()-h_ph))>3*TMath::Pi()/180){
                     continue;
                   }
                   //if(debug) cout << "passed phi" <<endl;                                                                
                   double diff = hypot(angle(mc.Phi()-pi2_ph)*sin(pi2_th),mc.Theta()-pi2_th);
                   //closest match which has negative charge                                                     
                   //if(debug) cout << "diff " << diff << endl;                                                            
                   if(diff < best_match_diff){
                     kbest2 = k;
                     best_match_diff = diff;
                   }
                 }
		 TLorentzVector pi1_truth, pi2_truth;

		 if(kbest1 >=0){
		   double pion_mass = db->GetParticle(211)->Mass();
		   pi1_truth = {mcparts->getPx(kbest1),mcparts->getPy(kbest1),mcparts->getPz(kbest1), 0};
		   pi1_truth.SetE(hypot(pi1_truth.P(),pion_mass));
		   pi1_truth_pid = mcparts->getPid(kbest1);
		   pi1_truth_z = pi1_truth.E()/(E-e_truth_p);
		   pi1_truth_p = pi1_truth.P();
		   pi1_truth_th = pi1_truth.Theta();
		   pi1_truth_ph = pi1_truth.Phi();
		   TLorentzVector pi1_truth_cm;
		   toCM(cm_truth, pi1_truth,pi1_truth_cm);
		   pi1_truth_cm_pt = pi1_truth_cm.Pt();
		   pi1_truth_cm_eta = pi1_truth_cm.Eta();
		   pi1_truth_cm_ph = pi1_truth_cm.Phi();
		 }
		 if(kbest2 >=0){
		   double pion_mass = db->GetParticle(211)->Mass();
                   pi2_truth_pid = mcparts->getPid(kbest2);
		   pi2_truth = {mcparts->getPx(kbest2),mcparts->getPy(kbest2),mcparts->getPz(kbest2), 0};
                   pi2_truth.SetE(hypot(pi2_truth.P(),pion_mass));
                   pi2_truth_z = pi2_truth.E()/(E-e_truth_p);
		   pi2_truth_p = pi2_truth.P();
		   pi2_truth_th = pi2_truth.Theta();
		   pi2_truth_ph = pi2_truth.Phi();
		   TLorentzVector pi2_truth_cm;
		   toCM(cm_truth, pi2_truth, pi2_truth_cm);
		   pi2_truth_cm_pt = pi2_truth_cm.Pt();
		   pi2_truth_cm_eta = pi2_truth_cm.Eta();
		   pi2_truth_cm_ph = pi2_truth_cm.Phi();
		 }
	       
		 if(dipion_tree != NULL){
		   if(kbest1>=0 && kbest2 >= 0){
		     diff_phi_cm_truth = angle(pi1_truth_cm_ph-pi2_truth_cm_ph);
		     diff_eta_cm_truth = pi1_truth_cm_eta-pi2_truth_cm_eta;
		     pair_mass_truth = (pi1_truth+pi2_truth).M();
		     mx_epi1pi2x_truth = (beam+target-el-pi1_truth-pi2_truth).M();
		     mx_epi1x_truth = (beam+target-el-pi1_truth).M();
		     mx_epi2x_truth = (beam+target-el-pi2_truth).M();
		   }
		   dipion_tree->Fill();
		 }
	       }
	     }
	   }
	 }

	 }
       }if(debug) cout << "end electron loop"<<endl;

       TLorentzVector mc;
       if(isMC){
	 for(int k = 0; k<mcparts->getRows();k++){
	   bool alreadymatched=0;
	   for(int j : matchedMCindices){
	     if(k == j)
	       alreadymatched=1;
	     //if(debug) cout << "skipping matched mc particle" << endl;
	   }
	   if(alreadymatched)
	     continue;
	   int pid = mcparts->getPid(k);
	   mc.SetXYZM(mcparts->getPx(k),mcparts->getPy(k),mcparts->getPz(k),db->GetParticle(pid)->Mass());
	   if(abs(pid) == 211 || abs(pid) == 321 || abs(pid) == 2212){

	     // only store unmatched generated hadrons in events where
	     // the electron is reconstructed
	     if(e_p == defval)
	       continue;
	     h_chi2pid = h_pid = h_p = h_th = h_ph = z = h_DC1x = h_DC2x = h_DC3x = h_DC1y = h_DC2y = h_DC3y = dvz = h_cm_th = h_cm_ph = h_cm_eta = h_cm_pt = h_cm_p = defval; 
	     h_truth_p = mc.P();
	     h_truth_th = mc.Theta();
	     h_truth_ph = mc.Phi();
	     h_truth_pid = pid;
	     if(hadron_tree != NULL)
	       hadron_tree->Fill();
	   }
	   
	   if(pid == 11 && -(mc-beam)*(mc-beam)> cut_Q2min && (target+beam-mc).M()> cut_Wmin && (E-mc.E())/E < cut_ymax){
	     e_p = e_th = e_ph = Q2 = nu = x = y = W = e_DC1x = e_DC2x = e_DC3x = e_DC1y = e_DC2y = e_DC3y = e_PCALx = e_PCALy = e_ecalfrac = e_pcal = e_vz = defval; 
	     e_truth_p = mc.P();
	     e_truth_th = mc.Theta();
	     e_truth_ph = mc.Phi();
	     e_truth_pid = 11;
	     
	     e_truth_vx = mcparts->getVx(k);
	     e_truth_vy = mcparts->getVy(k);
	     e_truth_vz = mcparts->getVz(k);
	     if(electron_tree != NULL)
	       electron_tree->Fill();
	   }
	 }
	 
       }
       
     } 
      
   }
   gBenchmark->Stop("timer");
   gBenchmark->Print("timer");
   //hm2gCut->SetLineColor(2);
   //hm2gCut->DrawCopy("same");
   TFile *f = new TFile(outputFile,"RECREATE");
   
   if(electron_tree != NULL)
     electron_tree->Write();
   if(hadron_tree != NULL)
     hadron_tree->Write();
   if(dihadron_tree != NULL)
     dihadron_tree->Write();
   if(dipion_tree != NULL)
     dipion_tree->Write();
   f->Close();
   
   
   

   auto finish = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = finish - start;
   std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<count<< " s\n";

}
