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
#include "DCfiducialcuts.h"

using namespace clas12;


const int INBEND=1, OUTBEND=0;

int debug = 0;
double defval = 0;

void SetLorentzVector(TLorentzVector &p4, clas12::region_part_ptr rp, double mass){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),mass);

}



double s30 = sin(TMath::Pi()/6);
double c30 = cos(TMath::Pi()/6);
/*bool dcOK(clas12::region_part_ptr p){
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
 }*/

bool pcalOK(clas12::region_part_ptr p){
   
  //cout << p->cal(PCAL)->getLv() << " " <<  p->cal(PCAL)->getLw() << " " <<  p->cal(PCAL)->getDv() << " " <<  p->cal(PCAL)->getDw() << " " <<  endl;
  return p->cal(PCAL)->getLv()>9.0 && p->cal(PCAL)->getLw()>9.0;
  /*double x = p->cal(PCAL)->getX();
  double y = p->cal(PCAL)->getY();

  double cutPCAL = 25; //in between sectors
  if (abs(x)<cutPCAL || abs(c30*y+s30*x)<cutPCAL || abs(c30*y-s30*x)<cutPCAL)
    return false;
  //cout << "pcal ok" <<endl;
  return true;*/
}

//list of cuts
double cut_ECfracmin = 0.17;
double cut_pcalmin = 0.07;
double cut_evzmin_inb =  -13; // inbending    -8;//-13;
double cut_evzmax_inb =  12;  // inbending     1;//12;

double cut_evzmin_outb = -18;
double cut_evzmax_outb = 10;

double cut_evzmin = cut_evzmin_inb;
double cut_evzmax = cut_evzmin_inb;


//double cut_Wmin  = 2;
double cut_Q2min = 1;
double cut_Wmin = 2;
double cut_ymax  = 0.85;

double cut_dvzmin = -20; //-5;//-20;
double cut_dvzmax = 20; //5;//20;

double cut_dtime_corr = 0.3;

//double cut_chi2pidmin = -2.5, cut_chi2pidmax = 2.5;

//maximum chi2pid difference for non-pions.  
//double cut_chi2pidmax_other = 3.0;


double cut_HTCCmin = 2;

double sf1[7], sf2[7], sf3[7], sf4[7];
double sfs1[7], sfs2[7], sfs3[7], sfs4[7];


//double cut_zmin = 0.3;

void toCM(TLorentzVector cm, TLorentzVector p,TLorentzVector& result){
  result = p;
  //  result.Print();
  result.RotateZ(TMath::Pi()-cm.Phi());
  result.RotateY(cm.Theta());

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
  int torus;
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
   TString qadbPath="";
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
    } else if (opt.Contains("--qadbPath=")){
      qadbPath=opt(11,opt.Sizeof());
      cout << "qadb path: "<< qadbPath <<endl;
    }

   }

   if(!(createDihadronTree||createHadronTree||createDipionTree||createElectronTree)){
     cout << "no trees selected!  Aborting" << endl;
     exit(0);
   }

   clas12::clas12databases *dbc12;
   if(!isMC){
     clas12::clas12databases::SetRCDBRemoteConnection();
     clas12::clas12databases::SetCCDBRemoteConnection();
     if(!qadbPath.EqualTo("")){
       cout << "setting qadb connection" <<endl;
       clas12::clas12databases::SetQADBConnection((const std::string)qadbPath);
       cout << "successfully connected to qadb" <<endl;
     }
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
   leaf(e_px);
   leaf(e_py);
   leaf(e_pz);
   leaf(nelectrons);
   leaf(nu);leaf(Q2);leaf(x);leaf(y);leaf(W);
   leaf(e_DC1x);leaf(e_DC2x);leaf(e_DC3x);leaf(e_DC1y);leaf(e_DC2y);leaf(e_DC3y);
   leaf(e_PCALx);leaf(e_PCALy);
   leaf(e_ecalfrac);leaf(e_pcal); leaf(e_ecalin); leaf(e_ecalout);
   leaf(e_vz);
   leaf(npip);leaf(npim);
   leaf(npp);leaf(npm);
   leaf(nKp);leaf(nKm);
   leaf(nh);
   leaf(ntracks);
   leaf(nhtracks);
   leaf(z_tot);

   
   // don't add these variables to the tree unless running MC 
   if(!isMC)
     tree = NULL;
   leaf(e_truth_pid);leaf(e_truth_p);leaf(e_truth_th);leaf(e_truth_ph);leaf(e_truth_px);leaf(e_truth_py);leaf(e_truth_pz);leaf(e_truth_vx);leaf(e_truth_vy); leaf(e_truth_vz);leaf(nhtracks_truth);
   

   //cout << "made leaves" <<endl;
   
   
   TTree* hadron_tree = createHadronTree ? new TTree("hadrons","hadrons") : NULL;
   tree = hadron_tree;
   leafx(nelectrons);
   leafx(E);leafx(helicity);leafx(e_p);leafx(e_th);leafx(e_ph);leafx(e_px);leafx(e_py);leafx(e_pz);
   leafx(nu);leafx(Q2);leafx(x);leafx(y);leafx(W);leafx(ntracks);leafx(nhtracks);
   leaf(h_chi2pid);leaf(h_pid);leaf(h_p);leaf(h_th);leaf(h_ph);leaf(h_px);leaf(h_py);leaf(h_pz);leaf(h_DC1x);leaf(h_DC1y);leaf(h_DC2x);leaf(h_DC2y);leaf(h_DC3x);leaf(h_DC3y);leaf(dvz);leaf(z); leaf(h_cm_p);leaf(h_cm_th);leaf(h_cm_ph);leaf(h_cm_eta);leaf(h_cm_rap);leaf(h_cm_pt);leaf(h_cm_zeta); leaf(h_eta); leaf(dtime); leaf(dtime_corr); leaf(missing_mass);

   
   if(!isMC) tree = NULL;
   leafx(e_truth_pid);leafx(e_truth_p);leafx(e_truth_th);leafx(e_truth_ph);leafx(e_truth_px);leafx(e_truth_py);leafx(e_truth_pz);
   leaf(h_truth_pid);leaf(h_truth_p);leaf(h_truth_th);leaf(h_truth_ph);leaf(h_truth_px);leaf(h_truth_py);leaf(h_truth_pz);leaf(h_truth_cm_p);leaf(h_truth_cm_th);leaf(h_truth_cm_ph);leaf(h_truth_cm_eta);leaf(h_truth_cm_rap);leaf(h_truth_cm_pt);leaf(h_truth_z);leaf(missing_mass_truth);leafx(nhtracks_truth);
   

   TTree* dihadron_tree = createDihadronTree ? new TTree("dihadrons","dihadrons") : NULL;
   tree = dihadron_tree;   
   leafx(E);leafx(helicity);leafx(e_p);leafx(e_th);leafx(e_ph);leafx(e_px);leafx(e_py);leafx(e_pz);leafx(nu);leafx(Q2);leafx(x);leafx(y);leafx(W); leafx(ntracks); leafx(nhtracks); leaf(diff_phi_cm_mix);leaf(diff_eta_cm_mix);leaf(diff_rap_cm_mix);leaf(h2_assoc_index);

leaf(pair_pt_cm);leaf(pair_phi_cm); leaf(pair_pt);leaf(pair_phi);
 

   // macro creates fields for two hadrons
#define leaf2(name) double h1_##name=0; double h2_##name=0; if(tree != NULL) {tree->Branch((TString)"h1_"+#name,&h1_##name,(TString)"h1_"+#name+(TString)"/D"); tree->Branch((TString)"h2_"+#name,&h2_##name,(TString)"h2_"+#name+(TString)"/D");}
   leaf2(chi2pid);leaf2(pid);leaf2(p);leaf2(th);leaf2(ph);leaf2(px);leaf2(py);leaf2(pz);leaf2(z);leaf2(eta);
   leaf2(cm_p);leaf2(cm_th);leaf2(cm_ph);leaf2(cm_eta);leaf2(cm_rap);leaf2(cm_pt);
   leaf(pair_mass);leaf(mx_eh1h2x);leaf(mx_eh1x);leaf(mx_eh2x);
   leaf(diff_eta);
   leaf(diff_phi);
   leafx(nelectrons);
   leaf(diff_phi_cm);
   leaf(diff_eta_cm);
   leaf(diff_rap_cm);
   leaf2(cm_zeta);

   
   if(!isMC) tree = NULL;
   leafx(e_truth_pid);leafx(e_truth_p);leafx(e_truth_th);leafx(e_truth_ph);leafx(e_truth_px);leafx(e_truth_py);leafx(e_truth_pz);
   leaf(h1_truth_z);
   leaf(h1_truth_pid);leaf(h1_truth_p);leaf(h1_truth_th);leaf(h1_truth_ph);leaf(h1_truth_px);leaf(h1_truth_py);leaf(h1_truth_pz);leaf(h1_truth_cm_p);
   leaf(h1_truth_cm_th);leaf(h1_truth_cm_ph);leaf(h1_truth_cm_eta);leaf(h1_truth_cm_rap);leaf(h1_truth_cm_pt);
   leaf(h2_truth_pid);leaf(h2_truth_p);leaf(h2_truth_th);leaf(h2_truth_ph);leaf(h2_truth_px);leaf(h2_truth_py);leaf(h2_truth_pz);leaf(h2_truth_cm_p);leaf(h2_truth_z);
   leaf(h2_truth_cm_th);leaf(h2_truth_cm_ph);leaf(h2_truth_cm_eta);leaf(h2_truth_cm_rap);leaf(h2_truth_cm_pt);
   leaf(diff_phi_cm_truth);leaf(diff_eta_cm_truth);leaf(diff_rap_cm_truth);leaf(pair_mass_truth);leaf(mx_eh1h2x_truth);leaf(mx_eh1x_truth);leaf(mx_eh2x_truth);leafx(nhtracks_truth);
   

   // A small tree for storing dipion events, without requiring one of them to be leading
   //
     
   TTree* dipion_tree = createDipionTree ? new TTree("dipions","dipions") : NULL;
   tree = dipion_tree;
   leafx(nelectrons);
   leafx(E);leafx(helicity);leafx(e_p);leafx(e_th);leafx(e_ph);leafx(nu);leafx(Q2);leafx(x);leafx(y);leafx(W);leafx(ntracks);leafx(nhtracks);
   leafx(h1_p);leafx(h2_p);leafx(h1_th);leafx(h2_th);leafx(h1_ph);leafx(h2_ph);
   leafx(h1_cm_eta);leafx(h2_cm_eta);leafx(h1_cm_rap);leafx(h2_cm_rap);leafx(h1_cm_pt);leafx(h2_cm_pt);leafx(h1_z);leafx(h2_z);leafx(h1_pid);leafx(h2_pid);leafx(h1_cm_ph);leafx(h2_cm_ph);

   leafx(diff_phi_cm);leafx(diff_eta_cm);leafx(pair_mass);leafx(mx_eh1h2x);leafx(mx_eh1x);leafx(mx_eh2x);
   
   if(!isMC) tree = NULL;
   leafx(e_truth_pid);leafx(e_truth_p);leafx(e_truth_th);leafx(e_truth_ph);
   
   leafx(h1_truth_pid);leafx(h1_truth_p);leafx(h1_truth_th);leafx(h1_truth_ph);leafx(h1_truth_z);
   leafx(h1_truth_cm_ph);leafx(h1_truth_cm_eta);leafx(h1_truth_cm_rap);leafx(h1_truth_cm_pt);
   leafx(h2_truth_pid);leafx(h2_truth_p);leafx(h2_truth_th);leafx(h2_truth_ph);leafx(h2_truth_z);
   leafx(h2_truth_cm_ph);leafx(h2_truth_cm_eta);leafx(h2_truth_cm_rap);leafx(h2_truth_cm_pt);
   leafx(diff_phi_cm_truth);leafx(diff_eta_cm_truth);leafx(diff_rap_cm_truth);leafx(pair_mass_truth);leafx(mx_eh1h2x_truth);leafx(mx_eh1x_truth);leafx(mx_eh2x_truth);
   
   

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
       if(!qadbPath.EqualTo("")){
	 c12.db().qadb_requireGolden(true); 
	 c12.applyQA();
       }
       clas12::rcdb_reader *rcdb = c12.rcdb();

       auto& current = rcdb->current();
       E = current.beam_energy/1000;
       
       clas12::ccdb_reader *ccdb = c12.ccdb();

       double torus_curr = current.torus_current;
       double torus_scale = current.torus_scale;
       if(torus_scale <0){
	 torus=INBEND;
	 cut_evzmin = cut_evzmin_inb;
	 cut_evzmax = cut_evzmax_inb;
       } else {
	 torus=OUTBEND;
	 cut_evzmin = cut_evzmin_outb;
	 cut_evzmax = cut_evzmax_outb;
	 cout << cut_evzmin << " " << cut_evzmax <<  endl;
       }


       cout << "torus current, scale = " << torus_curr << " " << torus_scale << endl;
       TString targetName(rcdb->current().target);
       cout << "target is " << targetName << "." << endl;
       std::cout << "beam energy is " << E << std::endl;
       rcdb->close();
       const TableOfDoubles_t& electron_sf = ccdb->requestTableDoubles("/calibration/eb/electron_sf");
       
       
       for(int row = 0; row<=6; row++){
	 sf1[row] = electron_sf[row][3];
	 sf2[row] = electron_sf[row][4];
         sf3[row] = electron_sf[row][5];
         sf4[row] = electron_sf[row][6];
         sfs1[row] = electron_sf[row][7];
         sfs2[row] = electron_sf[row][8];
         sfs3[row] = electron_sf[row][9];
         sfs4[row] = electron_sf[row][10];
       }
       ccdb->close();
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

     // for event mixing:  the cm frame fourvector two events ago
     TLorentzVector cm_mix = {0,0,0,0};
     // for event mixing:  the cm frame fourvector one events ago 
     TLorentzVector cm_mix_tmp = {0,0,0,0};
     // for event mixing: the trigger hadron in the previous entry.
     TLorentzVector had_mix = {0,0,0,0};

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
       
       if(c12.helonline() != NULL){
	 //cout << "helicity is not null" << endl;
	 //helicity = c12.helonline()->getHelicity();
	 helicity = c12.event()->getHelicity();
	 //cout << "helicity is " << helicity << endl;
       }
       else{
	 //cout << "helicity is NULL"<<endl;
       }
       //if(debug) cout << "helicity" << endl;
       TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
       
       mcpar_ptr mcparts;
       if(isMC){
	 mcparts = c12.mcparts();
	 nhtracks_truth =0;
	 //skip the first two entries; these are the initial state particles
	 //the intermediate state particles are removed by the pid cuts
	 for(int kk = 2; kk<mcparts->getRows();kk++){
	   int pid = mcparts->getPid(kk);
	   if(debug) cout << pid << " ";
	   if(abs(pid)==211 || abs(pid) == 321 || abs(pid)==2212)
	     nhtracks_truth++;
	 }
	 if(debug) cout << endl;
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
	 int sector = electrons[i]->getSector();
	 
	 e_DC1x=electrons[i]->traj(DC,DC1)->getX();
	 e_DC1y=electrons[i]->traj(DC,DC1)->getY();
	 e_DC2x=electrons[i]->traj(DC,DC3)->getX();
	 e_DC2y=electrons[i]->traj(DC,DC3)->getY();
	 e_DC3x=electrons[i]->traj(DC,DC6)->getX();
	 e_DC3y=electrons[i]->traj(DC,DC6)->getY();
	 if(!dcOK(electrons[i],torus) || !pcalOK(electrons[i]))
	   continue;
	 //if(debug) cout << "dc and pcal ok"<<endl;
	 //if(!dcok)
	 //  continue;
	 //if(debug) cout << "electron" << endl;

	 //the electron mass is a myth.  I could set this to zero and nothing would change.
	 //if(debug) cout << "CHECK 1" <<endl;
	 SetLorentzVector(el,electrons[i], 0.000511); 
	 e_p = el.P();
	 e_px = el.Px();
	 e_py = el.Py();
	 e_pz = el.Pz();
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

	 //double sf1 = 0, sf2 = 0, sf3 = 0, sf4 = 0;
	 //double sfs1 =0, sfs2=0, sfs3 = 0, sfs4 = 0;

	 //cout << sector << " " << sf1[sector] << " " <<sf2[sector] << " " <<sf3[sector] <<endl; 
	 
	 //ecal fraction
	 double ecalfrac_mu = sf1[sector]*(sf2[sector]+sf3[sector]/ecal+sf4[sector]/pow(ecal,2));
	 double ecalfrac_sigma = sfs1[sector]*(sfs2[sector]+sfs3[sector]/ecal+sfs4[sector]/pow(ecal,2));
	 double nsig = 3.5;
	 if(e_ecalfrac>ecalfrac_mu+nsig*ecalfrac_sigma || e_ecalfrac<ecalfrac_mu-nsig*ecalfrac_sigma)
	   continue;

	 if(e_pcal<cut_pcalmin)
	   continue;
	 e_ecalin = electrons[i]->cal(ECIN)->getEnergy();
	 e_ecalout = electrons[i]->cal(ECOUT)->getEnergy();
	 //cout << e_ecalin << " " << e_ecalout << " " << e_p << " " << e_pcal << endl;  
	 //->cal(ECOUT)->getEnergy(); 
	 
	 //further cut to remove pions
	 if(e_p> 4.5 && e_ecalin/e_p < 0.2 - e_pcal/e_p)
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
	 
	 if(Q2<cut_Q2min) 
	   continue;
	 if(W<cut_Wmin)
	   continue;
	 if(y > cut_ymax)
           continue;
	 
	 cm = beam+target-el;
	 
	 TLorentzVector target_cm;
	 toCM(cm, target,target_cm);
	 //check that the cm outgoing electron has phi=0
	 //TLorentzVector tmp;
	 //toCM(cm, el,tmp);
	 //cout << tmp.Phi() << endl;
	 //toCM(cm, beam-el, tmp);
	 //cout << "  " << tmp.Theta() << endl;

	 auto parts=c12.getDetParticles();
	 ntracks=0;
	 nhtracks=0;
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
	   if(abs(pid)==211 || abs(pid)==321 || abs(pid) == 2212)
	     nhtracks++;
	 }
	 if(debug) cout << ntracks << " tracks" << endl;
	 //if(debug) cout << "parts" << endl;
	 
	 double nHTCC = electrons[i]->che(HTCC)->getNphe();
	 if(nHTCC<=cut_HTCCmin)
	   continue;


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
	   e_truth_px = 0;
	   e_truth_py = 0;
	   e_truth_pz = 0;
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
	     e_truth_px = e_truth.Px();
	     e_truth_py = e_truth.Py();
	     e_truth_pz = e_truth.Pz();
	     cm_truth = beam+target-e_truth;
	     e_truth_vx = mcparts->getVx(kbest);
	     e_truth_vy = mcparts->getVy(kbest);
	     e_truth_vz = mcparts->getVz(kbest);
	   }
	 }
	 if(debug) cout << " check 1"<< endl;

	 
	 bool found_leader = 0;
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
	     h_DC2x=h->traj(DC,DC3)->getX();
	     h_DC2y=h->traj(DC,DC3)->getY();
	     h_DC3x=h->traj(DC,DC6)->getX();
	     h_DC3y=h->traj(DC,DC6)->getY();
	     
	     dtime = electrons[i]->getTime()-h->getTime();
	     
	     //if(debug) cout << "dc"<<endl;
	     if(!dcOK(h,torus))
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
	     
	     h_px = had.Px();
	     h_py = had.Py();
	     h_pz = had.Pz();

	     //h_pid = h->par()->getPid();
	     h_chi2pid = h->par()->getChi2Pid();
	     //if(debug) cout << "chi2pid" << endl;
	     //if(abs(h_chi2pid) > cut_chi2pidmax)
	     //continue;
	     
	     //scaling constant for pid cut
	     double C = .88*(h_pid == 211)+.93*(h_pid == -211)+ 1.30*(abs(h_pid)==2212)+1.0*(abs(h_pid) != 211 && abs(h_pid)!=2212);
	     if(abs(h_chi2pid) >3*C)
	       continue;
	     
	     //tighter cut for pions
	     if(abs(h_pid) == 211 && h2_p>2.44 && h_chi2pid <-(0.00869+14.98587*exp(-h_p/1.18236)+1.81751*exp(-h_p/4.86394))*C)
	       continue;

	     dvz = electrons[i]->par()->getVz()-h->par()->getVz();
	     
	     if(dvz < cut_dvzmin || dvz > cut_dvzmax)
	       continue;
	     
	     z = had.E()/nu;
	     
	     TLorentzVector h_cm;
	     toCM(cm, had,h_cm);
	     h_cm_p = h_cm.P();
	     h_cm_th = h_cm.Theta();
	     h_cm_eta = h_cm.PseudoRapidity();
	     h_cm_rap = h_cm.Rapidity();
	     h_cm_ph = h_cm.Phi();
	     h_cm_pt = h_cm.Pt();
	     h_cm_zeta = h_cm.E()/target_cm.E();
	     //if(debug) cout << h_cm_p << " " << pi_cm_th << " " << pi_cm_ph << endl;
	     
	     TLorentzVector h_truth;
	     int kbest = -1;
	     if(isMC){
	       double best_match_diff = 9999;
	       h_truth_pid = 0;
	       h_truth_p = 0;
	       h_truth_th = 0;
	       h_truth_ph = 0;
	       h_truth_px = 0;
               h_truth_py = 0;
               h_truth_pz = 0;
	       h_truth_cm_p = 0;
	       h_truth_cm_th = 0;
	       h_truth_cm_eta = 0;
	       h_truth_cm_rap = 0;
	       h_truth_cm_ph = 0;
	       h_truth_cm_pt = 0;
	       //int kbest = -1;
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
		 h_truth_px = h_truth.Px();
		 h_truth_py = h_truth.Py();
		 h_truth_pz = h_truth.Pz();
		 h_truth_z = h_truth.E()/(E-e_truth_p);
		 
		 h_truth_cm_p = h_truth_cm.P();
		 h_truth_cm_th = h_truth_cm.Theta();
		 h_truth_cm_eta = h_truth_cm.PseudoRapidity();
		 h_truth_cm_rap = h_truth_cm.Rapidity();
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
	     if((z > 0.4 && abs(h_pid)==211 && dihadron_tree!= NULL) || (dipion_tree != NULL && abs(h_pid)==211)){
	       h1_pid = h_pid;
	       h1_p = h_p;
	       h1_th = h_th;
	       h1_eta = h_eta;
	       h1_px = h_px;
               h1_py = h_py;
               h1_pz = h_pz;

	       h1_ph = h_ph;
	       h1_chi2pid = h_chi2pid;
	       
	       h1_cm_p = h_cm_p;
	       h1_cm_th = h_cm_th;
	       h1_cm_ph =h_cm_ph;
	       h1_cm_eta = h_cm_eta;
	       h1_cm_rap = h_cm_rap;
	       h1_cm_pt = h_cm_pt;
	       h1_cm_zeta = h_cm_zeta;
	       h1_z = z;
	       h1 = had;
	       if(isMC){
		 h1_truth_pid = h_truth_pid;
		 h1_truth_p = h_truth_p;
		 h1_truth_th = h_truth_th;
		 h1_truth_ph = h_truth_ph;
		 h1_truth_z = h_truth_z;
		 
		 h1_truth_cm_p = h_truth_cm_p;
		 h1_truth_cm_th = h_truth_cm_th;
		 h1_truth_cm_eta = h_truth_cm_eta;
		 h1_truth_cm_rap = h_truth_cm_rap;
		 h1_truth_cm_ph = h_truth_cm_ph;
		 h1_truth_cm_pt = h_truth_cm_pt;
	       }
	       h2_assoc_index = 0;
	       for(int k = 0; k<parts.size();k++){
		 if(k == j)
		   continue;






		 if(debug) cout << "check1" << endl;
		 auto h2 = parts[k];
		 int pid = h2->getPid();
		 h2_pid = pid;
		 if(h2_pid != 211 && h2_pid != -211 && h2_pid != 2212 && h2_pid != -2212 && h2_pid != 321 && h2_pid != -321)
		   continue;
		 TLorentzVector had2;
		 TLorentzVector had2_cm;
		 if(db->GetParticle(pid) == NULL || db->GetParticle(pid)->Mass()==0)
		   continue;
		 SetLorentzVector(had2,h2, db->GetParticle(pid)->Mass());
		 
		 toCM(cm, had2,had2_cm);
		 
		 h2_z = had2.E()/nu;
		 if(h2_z>h1_z)
		   continue;
		 h2_p = had2.P();
		 h2_th = had2.Theta();
		 h2_ph = had2.Phi();
		 h2_px = had2.Px();
		 h2_py = had2.Py();
		 h2_pz = had2.Pz();
		 h2_cm_pt = had2_cm.Pt();
		 h2_cm_ph = had2_cm.Phi();
		 h2_cm_eta = had2_cm.Eta();
		 h2_cm_rap = had2_cm.Rapidity();
		 h2_cm_zeta = had2_cm.E()/target_cm.M();
		 h2_pid = pid;
		 diff_phi_cm = angle(h1_cm_ph-h2_cm_ph);
		 diff_eta_cm = h1_cm_eta-h2_cm_eta;
		 diff_rap_cm = h1_cm_rap-h2_cm_rap;
		 
		 TLorentzVector pair_cm = had2_cm+h_cm;
		 pair_pt_cm = pair_cm.Pt();
		 pair_phi_cm = pair_cm.Phi();
		 
		 TLorentzVector pair = had2 + had;
		 pair_pt = pair.Pt();
                 pair_phi = pair.Phi();
		 

		 TLorentzVector had_cm_mix;
		 toCM(cm_mix, had,had_cm_mix);
		 diff_phi_cm_mix = angle(had_cm_mix.Phi()-h2_cm_ph);
		 diff_eta_cm_mix = had_cm_mix.Eta()-h2_cm_eta;
		 diff_rap_cm_mix = had_cm_mix.Rapidity()-h2_cm_rap;

		 pair_mass = (had2+had).M();
		 mx_eh1h2x = (beam+target-el-had-had2).M();
		 mx_eh1x = (beam+target-el-had).M();
		 mx_eh2x = (beam+target-el-had2).M();

		 dtime = electrons[i]->getTime()-h2->getTime();
		 
                 if(!dcOK(h2,torus))
                   continue;
                 double mass = db->GetParticle(h2_pid)->Mass();
                 SetLorentzVector(had2,h2, mass);
                 double c = 29.9792458; //cm/ns
                 dtime_corr =dtime-electrons[i]->getPath()/c+h2->getPath()/(had2.Beta()*c);
                 if(abs(dtime_corr) > cut_dtime_corr)
                   continue;
                 h2_chi2pid = h2->par()->getChi2Pid();

		 //scaling factor for chi2pid cuts
		 double C = .88*(h2_pid == 211)+.93*(h2_pid == -211)+ 1.30*(abs(h2_pid)==2212)+1.0*(abs(h2_pid) != 211 && abs(h2_pid)!=2212);
		 if(abs(h2_chi2pid) >3*C)
		   continue;
		 //tighter cut for pions
		 if(abs(h2_pid) == 211 && h2_p>2.44 && h2_chi2pid <-(0.00869+14.98587*exp(-h2_p/1.18236)+1.81751*exp(-h2_p/4.86394))*C)
		   continue;
                 //if(useCuts && abs(h2_chi2pid) > cut_chi2pidmax)
                 //  continue;

                 dvz = electrons[i]->par()->getVz()-h2->par()->getVz();

                 if(useCuts && (dvz < cut_dvzmin || dvz > cut_dvzmax))
                   continue;

		 // if only the dipion tree is being filled, and the second hadron is not a pion, skip this hadron
		 if(dihadron_tree == NULL && dipion_tree != NULL && abs(h2_pid) != 211)
		   continue;

		 
		 if(isMC){
		   h2_truth_z=0;
		   h2_truth_pid=0;
		   h2_truth_p=0;
		   h2_truth_th=0;
		   h2_truth_ph=0;
		   
		   h2_truth_px=0;
		   h2_truth_py=0;
		   h2_truth_pz=0;
		   h2_truth_cm_pt=0;
		   h2_truth_cm_ph=0;
		   h2_truth_cm_eta=0;
		   h2_truth_cm_rap=0;
		   
		   //double h1_p = had.P();
		   double h2_th = had2.Theta();
		   double h2_ph = had2.Phi();
		   double best_match_diff = 9999;
		   int kbest2 = -1;
		   
		   for(int kk = 0; kk<mcparts->getRows();kk++){
		     if(kk == kbest)
		       continue;
		     //if(debug) cout << "hadron" <<endl;                                                                                                                                        
		     TVector3 mc(mcparts->getPx(kk),mcparts->getPy(kk),mcparts->getPz(kk));
		     if(abs(mc.Theta()-h_th)>1*TMath::Pi()/180){
		       continue;
		     }
		     if(abs(angle(mc.Phi()-h_ph))>3*TMath::Pi()/180){
		       continue;
		     }
		     //if(debug) cout << "passed phi" <<endl;                                                                                                                                    
		     double diff = hypot(angle(mc.Phi()-h2_ph)*sin(h2_th),mc.Theta()-h2_th);
		     //closest match which has negative charge                                                                                                                                   
		     //if(debug) cout << "diff " << diff << endl;                                                                                                                                
		     if(diff < best_match_diff){
		       kbest2 = kk;
		       best_match_diff = diff;
		     }
		   }
		   TLorentzVector h1_truth = h_truth;
		   TLorentzVector h2_truth;
		   
		   if(kbest2 >=0){
		     h2_truth_pid = mcparts->getPid(kbest2);
		     double mass = db->GetParticle(h2_truth_pid)->Mass();
		     h2_truth = {mcparts->getPx(kbest2),mcparts->getPy(kbest2),mcparts->getPz(kbest2), 0};
		     h2_truth.SetE(hypot(h2_truth.P(),mass));
		     h2_truth_z = h2_truth.E()/(E-e_truth_p);
		     h2_truth_p = h2_truth.P();
		     h2_truth_th = h2_truth.Theta();
		     h2_truth_ph = h2_truth.Phi();

		     h2_truth_px = h2_truth.Px();
		     h2_truth_py = h2_truth.Py();
		     h2_truth_pz = h2_truth.Pz();
		     TLorentzVector h2_truth_cm;
		     toCM(cm_truth, h2_truth,h2_truth_cm);
		     h2_truth_cm_pt = h2_truth_cm.Pt();
		     h2_truth_cm_eta = h2_truth_cm.Eta();
		     h2_truth_cm_rap = h2_truth_cm.Rapidity();
		     h2_truth_cm_ph = h2_truth_cm.Phi();
		   }
		   if(isMC){
		     diff_phi_cm_truth = h2_truth_cm_ph-h1_truth_cm_ph;
		     diff_eta_cm_truth = h2_truth_cm_eta-h1_truth_cm_eta;
		   }
		 }
		 
		 if(dihadron_tree!= NULL && h1_z>0.4)
		   dihadron_tree->Fill();
		 if(dipion_tree != NULL && abs(h2_pid)==211)
		   dipion_tree->Fill();
		 h2_assoc_index++;
	       }

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
	     had_mix = had;//store the current hadron as the "mix" hadron.  
	   }
	 
	 
	 }
	 if(electron_tree != NULL)
           electron_tree->Fill();
         electrons_passCuts++;
	 cm_mix = cm_mix_tmp;
	 cm_mix_tmp = cm;
	   
       }
       if(debug) cout << "end electron loop"<<endl;
     
       //if(electron_tree != NULL)
       //electron_tree->Fill();
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
	     h_chi2pid = h_pid = h_p = h_th = h_ph = h_px = h_py = h_pz = z = h_DC1x = h_DC2x = h_DC3x = h_DC1y = h_DC2y = h_DC3y = dvz = h_cm_th = h_cm_ph = h_cm_eta = h_cm_pt = h_cm_p = defval; 
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
	     e_truth_px = mc.Px();
             e_truth_py = mc.Py();
             e_truth_pz = mc.Pz();

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
