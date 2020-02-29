/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "DDPlanarDigiProcessor.h"

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <EVENT/MCParticle.h>

#include <UTIL/CellIDEncoder.h>
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>


#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "AIDA/AIDA.h"

#include "marlin/ProcessorEventSeeder.h"
#include "marlin/AIDAProcessor.h"
#include "marlin/Global.h"

#include "CLHEP/Vector/TwoVector.h"

#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <climits>
#include <cfloat>

using namespace lcio ;
using namespace marlin ;
using namespace std ;


DDPlanarDigiProcessor aDDPlanarDigiProcessor ;

DDPlanarDigiProcessor::DDPlanarDigiProcessor() : Processor("DDPlanarDigiProcessor") {
  
  // modify processor description
  _description = "DDPlanarDigiProcessor creates TrackerHits from SimTrackerHits, smearing them according to the input parameters."
    "The geoemtry of the surface is taken from the DDRec::Surface asscociated to the hit via the cellID" ;
  
  
  // register steering parameters: name, description, class-variable, default value
  
  FloatVec resUEx ;
  resUEx.push_back( 0.0040 ) ;
  
  registerProcessorParameter( "ResolutionU" ,
                              "resolution in direction of u - either one per layer or one for all layers "  ,
                              _resU ,
                              resUEx) ;
  
  FloatVec resVEx ;
  resVEx.push_back( 0.0040 ) ;

  registerProcessorParameter( "ResolutionV" , 
                              "resolution in direction of v - either one per layer or one for all layers " ,
                             _resV ,
                              resVEx );

  registerProcessorParameter( "IsStrip",
                              "whether hits are 1D strip hits",
                              _isStrip,
                              bool(false) );
  
  
  registerProcessorParameter( "SubDetectorName" , 
                             "Name of dub detector" ,
                             _subDetName ,
                              std::string("VXD") );
    
  // Input collections
  registerInputCollection( LCIO::SIMTRACKERHIT,
                          "SimTrackHitCollectionName" , 
                          "Name of the Input SimTrackerHit collection"  ,
                          _inColName ,
                          std::string("VXDCollection") ) ;
  
  
  // Output collections
  registerOutputCollection( LCIO::TRACKERHITPLANE,
                           "TrackerHitCollectionName" , 
                           "Name of the TrackerHit output collection"  ,
                           _outColName ,
                           std::string("VTXTrackerHits") ) ;
  
  registerOutputCollection(LCIO::LCRELATION,
                           "SimTrkHitRelCollection",
                           "Name of TrackerHit SimTrackHit relation collection",
                           _outRelColName,
                           std::string("VTXTrackerHitRelations"));
  
  registerProcessorParameter( "ForceHitsOntoSurface" , 
                              "Project hits onto the surface in case they are not yet on the surface (default: false)" ,
                              _forceHitsOntoSurface ,
                              bool(false) );

  registerProcessorParameter( "MinimumEnergyPerHit" ,
                              "Minimum Energy (in GeV!) to accept hits, other hits are ignored",
                              _minEnergy,
                              double(0.0) );

  
  // setup the list of supported detectors
  
  
}

enum {
  hu = 0,
  hv,
  hitE,
  diffu,
  diffv,
  hSize 
} ;

void DDPlanarDigiProcessor::init() { 
  
  // usually a good idea to
  printParameters() ;
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
  // initialize gsl random generator
  _rng = gsl_rng_alloc(gsl_rng_ranlxs2);


  _h.resize( hSize ) ;

  Global::EVENTSEEDER->registerProcessor(this);

  
  if( _resU.size() !=  _resV.size() ) {
    
    std::stringstream ss ;
    ss << name() << "::init() - Inconsistent number of resolutions given for U and V coordinate: " 
       << "ResolutionU  :" <<   _resU.size() << " != ResolutionV : " <<  _resV.size() ;

    throw EVENT::Exception( ss.str() ) ;
  }

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();


  //===========  get the surface map from the SurfaceManager ================

  dd4hep::rec::SurfaceManager& surfMan = *theDetector.extension<dd4hep::rec::SurfaceManager>() ;

  dd4hep::DetElement det = theDetector.detector( _subDetName ) ;

  _map = surfMan.map( det.name() ) ;

  if( ! _map ) {   
    std::stringstream err  ; err << " Could not find surface map for detector: " 
                                 << _subDetName << " in SurfaceManager " ;
    throw Exception( err.str() ) ;
  }

  streamlog_out( DEBUG3 ) << " DDPlanarDigiProcessor::init(): found " << _map->size() 
                          << " surfaces for detector:" <<  _subDetName << std::endl ;

  streamlog_out( MESSAGE ) << " *** DDPlanarDigiProcessor::init(): creating histograms" << std::endl ;

  AIDAProcessor::histogramFactory(this) ; //->createHistogram1D( "hMCPEnergy", "energy of the MCParticles", 100 ) ;

  _h[ hu ] = new TH1F( "hu" , "smearing u" , 50, -5. , +5. );
  _h[ hv ] = new TH1F( "hv" , "smearing v" , 50, -5. , +5. );

  _h[ diffu ] = new TH1F( "diffu" , "diff u" , 1000, -5. , +5. );
  _h[ diffv ] = new TH1F( "diffv" , "diff v" , 1000, -5. , +5. );

  _h[ hitE ] = new TH1F( "hitE" , "hitEnergy in keV" , 1000, 0 , 200 );
  
}


void DDPlanarDigiProcessor::processRunHeader( LCRunHeader* ) {
  ++_nRun ;
} 

void DDPlanarDigiProcessor::processEvent( LCEvent * evt ) { 

  gsl_rng_set( _rng, Global::EVENTSEEDER->getSeed(this) ) ;   
  streamlog_out( DEBUG4 ) << "seed set to " << Global::EVENTSEEDER->getSeed(this) << std::endl;
  



  LCCollection* STHcol = 0 ;
  try{
    STHcol = evt->getCollection( _inColName ) ;
  }
  catch(DataNotAvailableException &e){
    streamlog_out(DEBUG4) << "Collection " << _inColName.c_str() << " is unavailable in event " << _nEvt << std::endl;
  }
  
  if( STHcol != 0 ){    
    


    unsigned nCreatedHits=0;
    unsigned nDismissedHits=0;
    
    LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHITPLANE )  ;
    
    CellIDEncoder<TrackerHitPlaneImpl> cellid_encoder( lcio::LCTrackerCellID::encoding_string() , trkhitVec ) ;

    LCCollectionVec* relCol = new LCCollectionVec(LCIO::LCRELATION);
    // to store the weights
    LCFlagImpl lcFlag(0) ;
    lcFlag.setBit( LCIO::LCREL_WEIGHTED ) ;
    relCol->setFlag( lcFlag.getFlag()  ) ;
    
    CellIDDecoder<SimTrackerHit> cellid_decoder( STHcol) ;
    
    
    int nSimHits = STHcol->getNumberOfElements()  ;
    
    streamlog_out( DEBUG4 ) << " processing collection " << _inColName  << " with " <<  nSimHits  << " hits ... " << std::endl ;
    
    for(int i=0; i< nSimHits; ++i){
      


      SimTrackerHit* simTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;

      _h[hitE]->Fill( simTHit->getEDep() * 1e6 );

      if( simTHit->getEDep() < _minEnergy ) {
        streamlog_out( DEBUG ) << "Hit with insufficient energy " << simTHit->getEDep()*1e6 << " keV" << std::endl;
        continue;
      }
      
      const int cellID0 = simTHit->getCellID0() ;
  
      //***********************************************************
      // get the measurement surface for this hit using the CellID
      //***********************************************************
      
      dd4hep::rec::SurfaceMap::const_iterator sI = _map->find( cellID0 ) ;

      if( sI == _map->end() ){    

        std::cout<< " DDPlanarDigiProcessor::processEvent(): no surface found for cellID : " 
                 <<   cellid_decoder( simTHit ).valueString() <<std::endl;

        
        std::stringstream err ; err << " DDPlanarDigiProcessor::processEvent(): no surface found for cellID : " 
                                    <<   cellid_decoder( simTHit ).valueString()  ;
        throw Exception ( err.str() ) ;
      }



      const dd4hep::rec::ISurface* surf = sI->second ;


      int layer  = cellid_decoder( simTHit )["layer"];



      dd4hep::rec::Vector3D oldPos( simTHit->getPosition()[0], simTHit->getPosition()[1], simTHit->getPosition()[2] );
      
      dd4hep::rec::Vector3D newPos ;

     //************************************************************
      // Check if Hit is inside senstive 
      //************************************************************
      
      if ( ! surf->insideBounds( dd4hep::mm * oldPos ) ) {
        
        streamlog_out( DEBUG3 ) << "  hit at " << oldPos 
                                << " " << cellid_decoder( simTHit).valueString() 
                                << " is not on surface " 
                                << *surf  
                                << " distance: " << surf->distance(  dd4hep::mm * oldPos )
                                << std::endl;        

        
        
        
        if( _forceHitsOntoSurface ){
          
          dd4hep::rec::Vector2D lv = surf->globalToLocal( dd4hep::mm * oldPos  ) ;
          
          dd4hep::rec::Vector3D oldPosOnSurf = (1./dd4hep::mm) * surf->localToGlobal( lv ) ; 
          
          streamlog_out( DEBUG3 ) << " moved to " << oldPosOnSurf << " distance " << (oldPosOnSurf-oldPos).r()
                                  << std::endl;        
            
          oldPos = oldPosOnSurf ;

        } else {

          ++nDismissedHits;
        
          continue; 
        }
      }
      
      //**************************************************************************
      // Try to smear the hit but ensure the hit is inside the sensitive region
      //**************************************************************************
      
      dd4hep::rec::Vector3D u = surf->u() ;
      dd4hep::rec::Vector3D v = surf->v() ;
      

      // get local coordinates on surface
      dd4hep::rec::Vector2D lv = surf->globalToLocal( dd4hep::mm * oldPos  ) ;
      double uL = lv[0] / dd4hep::mm ;
      double vL = lv[1] / dd4hep::mm ;


      bool accept_hit = false ;
      unsigned  tries   =  0 ;              
      static const unsigned MaxTries = 10 ; 
      
      float resU = ( _resU.size() > 1 ?   _resU.at(  layer )     : _resU.at(0)   )  ;
      float resV = ( _resV.size() > 1 ?   _resV.at(  layer )     : _resV.at(0)   )  ; 


      while( tries < MaxTries ) {
        
        if( tries > 0 ) streamlog_out(DEBUG0) << "retry smearing for " <<  cellid_decoder( simTHit ).valueString() << " : retries " << tries << std::endl;
        
        double uSmear  = gsl_ran_gaussian( _rng, resU ) ;
        double vSmear  = gsl_ran_gaussian( _rng, resV ) ;

        
        // dd4hep::rec::Vector3D newPosTmp = oldPos +  uSmear * u ;  
        // if( ! _isStrip )  newPosTmp = newPosTmp +  vSmear * v ;  
        
        
        dd4hep::rec::Vector3D newPosTmp = 1./dd4hep::mm  * ( ! _isStrip  ? 
                                                            surf->localToGlobal( dd4hep::rec::Vector2D (  ( uL + uSmear ) * dd4hep::mm, ( vL + vSmear )  *dd4hep::mm ) )  :
                                                            surf->localToGlobal( dd4hep::rec::Vector2D (  ( uL + uSmear ) * dd4hep::mm,          0.                  ) ) 
                                                            ) ;

        streamlog_out( DEBUG1 ) << " hit at    : " << oldPos 
                                << " smeared to: " << newPosTmp
                                << " uL: " << uL 
                                << " vL: " << vL 
                                << " uSmear: " << uSmear
                                << " vSmear: " << vSmear
                                << std::endl ;


        if ( surf->insideBounds( dd4hep::mm * newPosTmp ) ) {    
          
          accept_hit = true ;
          newPos     = newPosTmp ;

          _h[hu]->Fill(  uSmear / resU ) ; 
          _h[hv]->Fill(  vSmear / resV ) ; 

          _h[diffu]->Fill( uSmear );
          _h[diffv]->Fill( vSmear );

          break;  

        } else { 
          
          streamlog_out( DEBUG1 ) << "  hit at " << newPosTmp 
                                  << " " << cellid_decoder( simTHit).valueString() 
                                  << " is not on surface " 
                                  << " distance: " << surf->distance( dd4hep::mm * newPosTmp ) 
                                  << std::endl;        
        }
        
        ++tries;
       }
      
      
      if( accept_hit == false ) {
        streamlog_out(DEBUG4) << "hit could not be smeared within ladder after " << MaxTries << "  tries: hit dropped"  << std::endl;
        ++nDismissedHits;
        continue; 
      } 
      
      //**************************************************************************
      // Store hit variables to TrackerHitPlaneImpl
      //**************************************************************************
      

      TrackerHitPlaneImpl* trkHit = new TrackerHitPlaneImpl ;
                  
      const int cellID1 = simTHit->getCellID1() ;
      trkHit->setCellID0( cellID0 ) ;
      trkHit->setCellID1( cellID1 ) ;
      
      trkHit->setPosition( newPos.const_array()  ) ;
      trkHit->setTime( simTHit->getTime() ) ;
      trkHit->setEDep( simTHit->getEDep() ) ;

      float u_direction[2] ;
      u_direction[0] = u.theta();
      u_direction[1] = u.phi();
      
      float v_direction[2] ;
      v_direction[0] = v.theta();
      v_direction[1] = v.phi();
      
      streamlog_out(DEBUG0)  << " U[0] = "<< u_direction[0] << " U[1] = "<< u_direction[1] 
                             << " V[0] = "<< v_direction[0] << " V[1] = "<< v_direction[1]
                             << std::endl ;

      trkHit->setU( u_direction ) ;
      trkHit->setV( v_direction ) ;
      
      trkHit->setdU( resU ) ;

      if( _isStrip ) {

        // store the resolution from the length of the wafer - in case a fitter might want to treat this as 2d hit ....
        double stripRes = (surf->length_along_v() / dd4hep::mm ) / std::sqrt( 12. ) ;
        trkHit->setdV( stripRes ); 

      } else {
        trkHit->setdV( resV ) ;
      }

      if( _isStrip ){
        trkHit->setType( UTIL::set_bit( trkHit->getType() ,  UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ) ) ;
      }

      //**************************************************************************
      // Set Relation to SimTrackerHit
      //**************************************************************************    
         
      LCRelationImpl* rel = new LCRelationImpl;

      rel->setFrom (trkHit);
      rel->setTo (simTHit);
      rel->setWeight( 1.0 );
      relCol->addElement(rel);

      
      //**************************************************************************
      // Add hit to collection
      //**************************************************************************    
      
      trkhitVec->addElement( trkHit ) ; 
      
      ++nCreatedHits;
      
      streamlog_out(DEBUG3) << "-------------------------------------------------------" << std::endl;
      
    }      
    
    
    //**************************************************************************
    // Add collection to event
    //**************************************************************************    
    
    evt->addCollection( trkhitVec , _outColName ) ;
    evt->addCollection( relCol , _outRelColName ) ;
    
    streamlog_out(DEBUG4) << "Created " << nCreatedHits << " hits, " << nDismissedHits << " hits  dismissed as not on sensitive element\n";
    
  }
  _nEvt ++ ;
}



void DDPlanarDigiProcessor::check( LCEvent* ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void DDPlanarDigiProcessor::end(){ 

  gsl_rng_free( _rng );
  
  streamlog_out(MESSAGE) << " end()  " << name() 
  << " processed " << _nEvt << " events in " << _nRun << " runs "
  << std::endl ;
  
}
