#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

struct pho_obj {
  unsigned int idx;
  double pt;
};

// Initialize branches _____________________________________________________//
struct gen_obj {
//  unsigned int idx;
  double pt;
  double eta;
//  double phi;
};
std::vector<gen_obj> vA_;

// Initialize branches _____________________________________________________//
void SCRegressor::branchesDiPhotonSel ( TTree* tree, edm::Service<TFileService> &fs )
{
  tree->Branch("m0",        &m0_);
  tree->Branch("FC_inputs", &vFC_inputs_);
  tree->Branch("hltAccept", &hltAccept_);
  tree->Branch("nRecoPho",  &nRecoPho_);
  tree->Branch("minDR",     &vMinDR_);
  tree->Branch("mA_close",     &mA_);
 // tree->Branch("mA_Recoclose",     &mA2_);
  tree->Branch("mA_closeGen",     &mA1_);
  tree->Branch("GenMatch",     &vGenMatch_);

 //hNpassed_nKin      = fs->make<TH1F>("hNpassed_nRecoPho", "isPassed;isPassed;N", 2, 0., 2);
  hNpassed_hlt      = fs->make<TH1F>("hNpassed_hlt", "isPassed;isPassed;N", 2, 0., 2);
  hNpassed_nRecoPho = fs->make<TH1F>("hNpassed_nRecoPho", "isPassed;isPassed;N", 2, 0., 2);
  hNpassed_presel   = fs->make<TH1F>("hNpassed_presel", "isPassed;isPassed;N", 2, 0., 2);
  hNpassed_mGG      = fs->make<TH1F>("hNpassed_mGG", "isPassed;isPassed;N", 2, 0., 2);
  hNpassed_pt0mGG   = fs->make<TH1F>("hNpassed_pt0mGG", "isPassed;isPassed;N", 2, 0., 2);
  hFill_PairMatch   = fs->make<TH1F>("hFill_PairMatch", "isPassed;isPassed;N", 2, 0., 2);
//  hFill_Gmatch      =fs->make<TH1F>("hFill_Gmatch", "isPassed;isPassed;N",2, 0., 2);
  hFill_Amatch      =fs->make<TH1F>("hFill_Amatch", "isPassed;isPassed;N",2, 0., 2);
  hFill_Gmatch      =fs->make<TH1F>("hFill_Gmatch", "isPassed;isPassed;N",2, 0., 2);
}



// Run event selection ___________________________________________________________________//
bool SCRegressor::runDiPhotonSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  //edm::Handle<PhotonCollection> photons;
 // iEvent.getByToken(photonCollectionT_, photons);


  edm::Handle<edm::TriggerResults> trgs;
  iEvent.getByToken( trgResultsT_, trgs );
  //std::cout<<std::endl;
  const edm::TriggerNames &triggerNames = iEvent.triggerNames( *trgs );
  if (debug)  std::cout << " N triggers:" << trgs->size() << std::endl;
  for ( unsigned int iT = 0; iT < trgs->size(); iT++ ) {
   if (debug)  std::cout << " name["<<iT<<"]:"<<triggerNames.triggerName(iT)<< std::endl;
  }

  int hltAccept = -1;
  std::string trgName = "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_*_Mass55_v*";
//  std::string trgName = "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_*Mass90_v*";
  std::vector< std::vector<std::string>::const_iterator > trgMatches = edm::regexMatch( triggerNames.triggerNames(), trgName );
  if (debug) std::cout << " N matches: " << trgMatches.size() << std::endl;
   if ( !trgMatches.empty() ) {
    hltAccept = 0;
    for ( auto const& iT : trgMatches ) {
      if ( debug ) std::cout << " name["<<triggerNames.triggerIndex(*iT)<<"]:"<< *iT << " -> " << trgs->accept(triggerNames.triggerIndex(*iT)) << std::endl;
    if ( trgs->accept(triggerNames.triggerIndex(*iT)) ) hltAccept = 1;
      break;
    }
  }
  hltAccept_ = hltAccept;
  hNpassed_hlt->Fill(0.);
  if ( hltAccept_ != 1 ) return false;
  hNpassed_hlt->Fill(1.);




  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);
  ////////// Apply selection on reco photons  //////////

  if ( debug ) std::cout << " Pho collection size:" << photons->size() << std::endl;

   hNpassed_nRecoPho->Fill(0.);
  // Count number of "reco" photons
  std::vector<unsigned int> vRecoPhoIdxs;
  for ( unsigned int iP = 0; iP < photons->size(); iP++ ) {
    PhotonRef iPho( photons, iP );
    //if ( std::abs(iPho->pt()) < 5. ) continue;
    if ( std::abs(iPho->pt()) < 10. ) continue;
    vRecoPhoIdxs.push_back( iP );
  }
  if ( debug ) std::cout << " Reco pho size:" << vRecoPhoIdxs.size() << std::endl;
  nRecoPho_ = vRecoPhoIdxs.size();

  // Ensure at least 2 kinematic trigger-like photons
  /*std::vector<unsigned int> vKinPhoIdxs;
  for ( unsigned int iP : vRecoPhoIdxs ) {
    PhotonRef iPho( photons, iP );
    if ( std::abs(iPho->pt()) <= 18. ) continue;
    //if ( std::abs(iPho->pt()) < 20. ) continue;
    if ( std::abs(iPho->eta()) >= 1.442 ) continue;
    //if ( std::abs(iPho->eta()) >= 2.4 ) continue;
    vKinPhoIdxs.push_back( iP );
  }*/
  //if ( vRecoPhoIdxs.size() < 2 ) return false;

 // if ( vRecoPhoIdxs.size() != 2 && vRecoPhoIdxs.size() != 3 ) return false; // 2 or 3 for Mike.
  if ( vRecoPhoIdxs.size() != 3 ) return false;//exactly 3 for Abhi.
  hNpassed_nRecoPho->Fill(1.);

   /// Apply preselection criteria ////


  hNpassed_presel->Fill(0.);
  std::vector<pho_obj> vPhos;
  for ( unsigned int iP : vRecoPhoIdxs ) {

    PhotonRef iPho( photons, iP );

    ///*
    if ( iPho->full5x5_r9() <= 0.5 ) continue;
    if ( iPho->hadTowOverEm() >= 0.08 ) continue;
    if ( iPho->hasPixelSeed() == true ) continue;
    //if ( iPho->passElectronVeto() == true ) continue;
    //if ( iPho->userFloat("phoChargedIsolation")/std::abs(iPho->pt()) > 0.3 ) continue;
    //*/

    ///*
    if ( iPho->full5x5_r9() <= 0.85 ) {
      if ( iPho->full5x5_sigmaIetaIeta() >= 0.015 ) continue;
     // if ( iPho->userFloat("phoPhotonIsolation") >= 4.0 ) continue;
      if ( iPho->photonIso() >= 4.0 ) continue;
      if ( iPho->trkSumPtHollowConeDR03() >= 6. ) continue;
      //if ( iPho->trackIso() >= 6. ) continue;
    }
    //*/
    if ( debug ) std::cout << " >> pT:" << iPho->pt() << " eta:" << iPho->eta() << " phi: " << iPho->phi() << " E:" << iPho->energy() << std::endl;

    pho_obj Pho_obj = { iP, std::abs(iPho->pt()) };
    vPhos.push_back( Pho_obj );

  } // preselected photons
  if ( debug ) std::cout << " Presel pho size:" << vPhos.size() << std::endl;

//   if ( vPhos.size() != 2 ) return false; //Ensure 2 presel photon for Mike

   if ( vPhos.size() != 3 ) return false; //Ensure 3 presel photon for Abhi
  hNpassed_presel->Fill(1.);

  // Sort photons by pT, for abitrary N
  std::sort( vPhos.begin(), vPhos.end(), [](auto const &a, auto const &b) { return a.pt > b.pt; } );
  for ( unsigned int iP = 0; iP < vPhos.size(); iP++ ) {
    PhotonRef iPho( photons, vPhos[iP].idx );
    if ( debug ) std::cout << " >> pT:" << iPho->pt() << " eta:" << iPho->eta() << " phi: " << iPho->phi() << " E:" << iPho->energy() << std::endl;
  }

  // Check if any photon pairing passes invariant mass cut
  hNpassed_mGG->Fill(0.);
  std::vector<int> vPhoIdxs;
  bool passedMassCut = false;

	//Mike///
/* for ( unsigned int j = 0; j < vPhos.size()-1; j++ ) {

    PhotonRef jPho( photons, vPhos[j].idx );

    for ( unsigned int k = 1; k < vPhos.size(); k++ ) {

      if ( k <= j ) continue;
      PhotonRef kPho( photons, vPhos[k].idx );
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vDiPho = jPho->p4() + kPho->p4();
      if ( debug ) std::cout << " >> m0:" << vDiPho.mass() << std::endl;

      //if ( vDiPho.mass() >= 100. && vDiPho.mass() <= 180.) {
      if ( vDiPho.mass() > 90. ) {
        vPhoIdxs.push_back( vPhos[j].idx );
        vPhoIdxs.push_back( vPhos[k].idx );
        m0_ = vDiPho.mass();
        passedMassCut = true;
        break;
       }
      } //k
     if ( passedMassCut ) break;
    } //j
  */
//////// Abhirami //////
   PhotonRef jPho( photons, vPhos[0].idx );
     PhotonRef kPho( photons, vPhos[1].idx );
     PhotonRef lPho( photons, vPhos[2].idx );

    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vTriPho = jPho->p4() + kPho->p4() + lPho->p4();

	if (vTriPho.mass() > 90)   {
        vPhoIdxs.push_back( vPhos[0].idx );
        vPhoIdxs.push_back( vPhos[1].idx );
        vPhoIdxs.push_back( vPhos[2].idx );
        m0_ = vTriPho.mass();
	passedMassCut = true;

     }

  if ( !passedMassCut ) return false;
  if ( debug ) std::cout << " >> m0:" << m0_ << std::endl;
  hNpassed_mGG->Fill(1.);

  // Apply diphoton pT cuts
  //
  hNpassed_pt0mGG->Fill(0.);

//  float ptCut[2]   = { 30., 18.}; //MIke
 // float ptOmCut[2] = {  3.,  4. };

  float ptCut[3]   = { 30., 18. , 0.}; //Abhi
  float ptOmCut[3] = {  3.,  4. , 100.};
  vPreselPhoIdxs_.clear();


   for ( unsigned int iP = 0; iP < vPhoIdxs.size(); iP++ ) {

    PhotonRef iPho( photons, vPhoIdxs[iP] );

    if ( std::abs(iPho->pt()) < ptCut[iP] ) continue;
    if ( std::abs(iPho->pt()) < m0_/ptOmCut[iP] ) continue;
    if ( debug ) std::cout << " >> pT:" << iPho->pt() << " eta:" << iPho->eta() << " phi: " << iPho->phi() << " E:" << iPho->energy() << std::endl;

    vPreselPhoIdxs_.push_back( vPhoIdxs[iP] );

  } // vPhoIdxs


  // ****  Finding closest dR *********
  float pdR, mdR;

  unsigned int id0, id1;
  mdR=50.;

  id0=0, id1=0; //id2=10;

  for (unsigned int iP = 0; iP < vPhos.size()-1; iP++ ) {

    PhotonRef iPho( photons, vPhos[iP].idx );

    if ( std::abs(iPho->pt()) < ptCut[iP] ) continue;
    if ( std::abs(iPho->pt()) < m0_/ptOmCut[iP] ) continue;
    if ( debug ) std::cout << " >> pT:" << iPho->pt() << " eta:" << iPho->eta() << " phi: " << iPho->phi() << " E:" << iPho->energy() << std::endl;

    for (unsigned int gP = 1; gP < vPhos.size(); gP++ ) {

      if (gP <= iP) continue;
      PhotonRef gPho( photons, vPhos[gP].idx );
      if ( std::abs(gPho->pt()) < ptCut[gP] ) continue;
      if ( std::abs(gPho->pt()) < m0_/ptOmCut[gP] ) continue;

      pdR = std::abs(reco::deltaR(iPho->eta(),iPho->phi(), gPho->eta(), gPho->phi()));
      //   std::cout<<"\n Indices: "<<iP<<" "<<gP<<" dR= "<<pdR<<std::endl;
      if (pdR < mdR) {

        mdR = pdR;
        id0 = vPhos[iP].idx;
        id1 = vPhos[gP].idx;
        std::cout<<"\n Entered Indices: "<<id0<<" "<<id1<<" dR= "<<pdR<<std::endl;
        /*             for(unsigned int lP =0; lP <vPhos.size(); lP++){
                      if (vPhos[lP].idx != id0 && vPhos[lP].idx != id1)
                                  id2= vPhos[lP].idx;

           }
        */
        //    vPreselPhoIdxs_.push_back( vPhos[iP].idx );

      } //dR id loop

    }//gP

  }//iP

  ///std::cout<<"\n Entered"<<std::endl;
  /// std::cout<<"\n id0, id1, id2= "<<id0<<" "<<id1<<" "<<id2<<std::endl;
  ///if( id2!=10)

  /* This is segfaulting
  unsigned int id2 = 3-id0-id1;
  ///std::cout<<"\n id0, id1, id2= "<<id0<<" "<<id1<<" "<<id2<<std::endl;
  if (id0!=id1)  vPreselPhoIdxs_.push_back( vPhos[id2].idx );
  */


  //  if ( vPreselPhoIdxs_.size() != 2 ) return false; //Mike
  //   ******  //
  //   if ( vPreselPhoIdxs_.size() != 3 ) return false; //Abhi

  if ( debug ) std::cout << " Reco pho size:" << vPhos.size() << std::endl;
  if ( debug ) std::cout << " >> Passed selection. " << std::endl;
  hNpassed_pt0mGG->Fill(1.);

  /*
  // Ensure exactly two "reco" photons
  hNpassed_nRecoPho->Fill(0.);
  if ( nRecoPho_ != 2 ) return false;
  hNpassed_nRecoPho->Fill(1.);
  */
  return true;
}
// Fill branches ___________________________________________________________________//
void SCRegressor::fillDiPhotonSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup )

{

 //std::cout<<"\nenter fillDi"<<std::endl;
 // count ++;
 hFill_PairMatch->Fill(0.);
  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

 std::vector<unsigned int> vGenMatchIdxs_;
 std::vector<unsigned int> vGenMatchIdxsAll_;
 vGenMatchIdxs_.clear();
 mA_.clear();
 mA1_.clear();
 mA2_.clear();
 vGenMatch_.clear();
//vA_.clear();
  float gdR, gdR_, recoDR, recoDR_;
  int recoDR_idx, recoDR_idx_;

  // Fill kinematic variables
  float dphi[2] = { 0., 0. };
  //float dphi[3] = { 0., 0., 0. };
  vFC_inputs_.clear();
  for ( unsigned int iP = 0; iP < vRegressPhoIdxs_.size(); iP++ ) {
    PhotonRef iPho( photons, vRegressPhoIdxs_[iP] );
    vFC_inputs_.push_back( iPho->pt()/m0_ );
    vFC_inputs_.push_back( iPho->eta() );
    dphi[iP] = iPho->phi();
  }
  vFC_inputs_.push_back( TMath::Cos(reco::deltaPhi(dphi[0], dphi[1])) );
}

/*

 //Find the closest preselected photons



 std::vector<unsigned int> vClosestdRidxs_;

  float pdR, mdR ;
 unsigned int id0, id1;//, id2;
 // id0=-1, id1=-1, id2=-1;
  vClosestdRidxs_.clear();
  vMinDR_.clear();
  mdR=50.;
  for ( unsigned int fP = 0; fP < vPreselPhoIdxs_.size()-1; fP++ ) {
      PhotonRef fPho( photons, vPreselPhoIdxs_[fP] );

    for ( unsigned int gP = 1; gP < vPreselPhoIdxs_.size(); gP++ ) {
       	if (gP <= fP) continue;
 	PhotonRef gPho( photons, vPreselPhoIdxs_[gP] );

	pdR = std::abs(reco::deltaR(fPho->eta(),fPho->phi(), gPho->eta(), gPho->phi()));
//	std::cout<<"\n Indices: "<<fP<<" "<<gP<<" dR= "<<pdR<<std::endl;
  	if (pdR < mdR) {
		mdR = pdR;
		id0 = vPreselPhoIdxs_[fP];
		id1 = vPreselPhoIdxs_[gP];
*/
        /*	for( unsigned int lP =0; lP <vPreselPhoIdxs_.size(); lP++){
		    while (vPreselPhoIdxs_[lP] != id0 && vPreselPhoIdxs_[lP] != id1)
				id2= vPreselPhoIdxs_[lP];
		}*/
/*	}
   }
 }
 vClosestdRidxs_.push_back(id0);
 vClosestdRidxs_.push_back(id1);
 //vClosestdRidxs_.push_back(id2);
 vMinDR_.push_back( mdR );
// std::cout<<"\nsize of vector= "<<vClosestdRidxs_.size()<<std::endl;
if (debug) std::cout<<"\n \n Closest indices: " <<std::endl;
if (debug) std::cout<<vClosestdRidxs_.at(0)<<" "<<vClosestdRidxs_.at(1)<<std::endl;//" "<<vClosestdRidxs_.at(2)<<std::endl;
*/
 // ***  Get Gen matching of preselected photons *** ///
/*
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

// not for bkg ////////

//std::cout<<"\n Gen part size: "<<genParticles->size()<<std::endl;

 vGenMatchIdxsAll_.clear();
for ( unsigned int p = 0; p < vPreselPhoIdxs_.size(); p++ ) {
     PhotonRef pPho (photons, vPreselPhoIdxs_[p]);

     recoDR_=0.04;
     recoDR_idx_=-1;
     for ( unsigned int iP = 0; iP < genParticles->size(); iP++ ) {

        reco::GenParticleRef pGen( genParticles, iP );
        if ( std::abs(pGen->pdgId()) != 22 ) continue;

        gdR_ = std::abs(reco::deltaR(pGen->eta(),pGen->phi(), pPho->eta(),pPho->phi()));
        if ( gdR_ < recoDR_ ) {
          recoDR_ = gdR_;
          recoDR_idx_ = iP;
     //   std::cout<<"\n recodr_idx at "<<iG<<" =" <<recoDR_idx<<std::endl;
         }
    }
  	 if (recoDR_idx_ == -1) continue;
        vGenMatchIdxsAll_.push_back(recoDR_idx_);
}
   //reco::GenParticleRef kpGen( genParticles, vPreselPhoIdxs_[p] );
   int Acount=0;
for (unsigned int p=0; p< vGenMatchIdxsAll_.size(); p++)
   {    reco::GenParticleRef pGen( genParticles, vGenMatchIdxsAll_[p] );
        hFill_Amatch->Fill(0.);
        if ( std::abs(pGen->mother()->pdgId()) != 25 ) continue;
        if (pGen->mother()->numberOfDaughters() != 2 ) continue;
//        gen_obj A_Obj = { std::abs(kGen->mother()->pt()),std::abs(kGen->mother()->eta())};
	Acount++;
        hFill_Amatch->Fill(1.);

   }
   hFill_Gmatch->Fill(0.);
   if (Acount ==1)
        { hFill_Gmatch->Fill(1.);
	  vGenMatch_.push_back(1.);}
     else
	vGenMatch_.push_back(0.);
   //if (Acount ==3) hFill_Gmatch->Fill(1.);
}
*/
 //  std::cout<<"\n\nSize of Acount: "<<Acount<<std::endl;

// not for bkg

// not for bkg
/*
////  GenMatching of closest presel //////
PhotonRef Pho1( photons, vClosestdRidxs_[0] );
PhotonRef Pho2( photons, vClosestdRidxs_[1] );
ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vClosPho = Pho1->p4() + Pho2->p4();
//std::cout << "\n\n >> Photon p4s" << Pho1->p4() <<" "<<Pho2->p4() << std::endl;
//std::cout << " >> close mA:" << vClosPho.mass() << std::endl;
mA_.push_back(vClosPho.mass());
*/

/*
////////////////
//std::cout<<"\n Gen particles"<<std::endl;
for ( unsigned int k = 0; k < vClosestdRidxs_.size(); k++ ) {
     PhotonRef kPho (photons, vClosestdRidxs_[k]);

     recoDR=0.04;
     recoDR_idx=-1;
     for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {

    	reco::GenParticleRef iGen( genParticles, iG );
	if ( std::abs(iGen->pdgId()) != 22 ) continue;
	//std::cout<<"\n Gen particles"<<std::endl;
//	std::cout<<"\n eta: "<<iGen->eta()<<" p4: "<<iGen->p4()<<std::endl;
        gdR = std::abs(reco::deltaR(iGen->eta(),iGen->phi(), kPho->eta(),kPho->phi()));
      	if ( gdR < recoDR ) {
          recoDR = gdR;
          recoDR_idx = iG;
  //      std::cout<<"\n recodr=  "<<recoDR<<" recodr_idx at "<<iG<<" =" <<recoDR_idx<<std::endl;
	 }
	}
//   hFill_Gmatch->Fill(0.);
 //  std::cout<<"\n gdR = "<<gdR<<" reco dR= "<<recoDR<<std::endl;
 //  std::cout<<"\n recodr_idx  =" <<recoDR_idx<<std::endl;
 if (recoDR_idx == -1) continue;
//   hFill_Gmatch->Fill(1.);
//   std::cout<<"\n gdR = "<<gdR<<" reco dR= "<<recoDR<<std::endl;
//   std::cout<<"\n recodr_idx  =" <<recoDR_idx<<std::endl;
   if ( std::find(vGenMatchIdxs_.begin(), vGenMatchIdxs_.end(), recoDR_idx) != vGenMatchIdxs_.end() ) continue;
//      vGenMatchIdxs_.push_back( minDR_idx );
   vGenMatchIdxs_.push_back (recoDR_idx) ;

}

//std::cout<<"\ngen matched"<<std::endl;
//std::cout<<"\ngen size "<<vGenMatchIdxs_.size()<<std::endl;
//std::cout<<vGenMatchIdxs_.at(0)<<" "<<vGenMatchIdxs_.at(1)<<std::endl;
  //hFill_PairMatch->Fill(1.);

//////  Not for bkg  /////
*/
/*
vA_.clear();

for (unsigned int iK =0; iK < vGenMatchIdxs_.size(); iK++){

	reco::GenParticleRef kGen( genParticles, vGenMatchIdxs_[iK] );
 	//const reco::Candidate* kGen->mother() = kGen->mother();
 //	hFill_Amatch->Fill(0.);
        if ( std::abs(kGen->mother()->pdgId()) != 25 ) continue;
        if (kGen->mother()->numberOfDaughters() != 2 ) continue;
	gen_obj A_Obj = { std::abs(kGen->mother()->pt()),std::abs(kGen->mother()->eta())};//, std::abs(kGen->mother()->phi()) };
        vA_.push_back(A_Obj);
//	hFill_Amatch->Fill(1.);

  //      std::cout<<"\nvA_ size " <<vA_.size()<<std::endl;
	 //if( vA_.size()==0 ) continue;
        //std::cout<<vA_.at(0).idx<<std::endl;

}
  //std::cout<<"\n Enter matching"<<std::endl;
 // hFill_PairMatch->Fill(0.);
  //std::cout<<"\nvA_ size " <<vA_.size()<<std::endl;
  //hFill_PairMatch->Fill(2.);
  if( vA_.size()==2 )

//  {std::cout<<vA_.at(0)<<" "<<vA_.at(1)<<" \n";}
 // {std::cout<<vA_.at(0).pt<<" "<<vA_.at(1).pt<<" \n";}
    { */


    //	hFill_PairMatch->Fill(1.);


//// not for bkg ///////
/*
 if (vGenMatchIdxs_.size()==2)
  {
	 reco::GenParticleRef Gen1( genParticles, vGenMatchIdxs_[0] );
       	 reco::GenParticleRef Gen2( genParticles, vGenMatchIdxs_[1] );
         ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vClosGen = Gen1->p4() + Gen2->p4();
        //std::cout<<"\n\n Gen1 p4: "<<Gen1->p4()<<" Gen 2 p4: "<<Gen2->p4()<<std::endl;
	//std::cout<<"Sum p4: "<<vClosGen<<" mass: "<<vClosGen.mass()<<std::endl;
	mA1_.push_back(vClosGen.mass());
  }
//      mA2_.push_back(vClosPho.mass());
*/
/*	if (vA_.at(0).pt==vA_.at(1).pt && vA_.at(0).eta==vA_.at(1).eta)// && vA_.at(0).phi==vA_.at(1).phi)
		{ hFill_PairMatch->Fill(2.);}
    }
*/

//}

