#include "AllpixDigiProcessor.h"

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


AllpixDigiProcessor aAllpixDigiProcessor ;

AllpixDigiProcessor::AllpixDigiProcessor() : Processor("AllpixDigiProcessor") {
  
  // modify processor description
  _description = "AllpixDigiProcessor creates TrackerHits from SimTrackerHits, smearing them according to the input parameters."
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

void AllpixDigiProcessor::init() {
  
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

  streamlog_out( DEBUG3 ) << " AllpixDigiProcessor::init(): found " << _map->size()
                          << " surfaces for detector:" <<  _subDetName << std::endl ;

  streamlog_out( MESSAGE ) << " *** AllpixDigiProcessor::init(): creating histograms" << std::endl ;

  AIDAProcessor::histogramFactory(this) ; //->createHistogram1D( "hMCPEnergy", "energy of the MCParticles", 100 ) ;

  _h[ hu ] = new TH1F( "hu" , "smearing u" , 50, -5. , +5. );
  _h[ hv ] = new TH1F( "hv" , "smearing v" , 50, -5. , +5. );

  _h[ diffu ] = new TH1F( "diffu" , "diff u" , 1000, -5. , +5. );
  _h[ diffv ] = new TH1F( "diffv" , "diff v" , 1000, -5. , +5. );

  _h[ hitE ] = new TH1F( "hitE" , "hitEnergy in keV" , 1000, 0 , 200 );
  
}


void AllpixDigiProcessor::processRunHeader( LCRunHeader* ) {
  ++_nRun ;
} 

void AllpixDigiProcessor::processEvent( LCEvent * evt ) {

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


    int nSimHits = STHcol->getNumberOfElements()  ;

    std::cout << " processing collection " << _inColName  << " with " <<  nSimHits  << " hits ... " << std::endl ;

    for(int i=0; i< nSimHits; ++i){

      SimTrackerHit* simTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;

      std::cout << simTHit->getPosition()[0]<<"\t"<< simTHit->getPosition()[1]<<"\t"<< simTHit->getPosition()[2] <<std::endl;

      }

  _nEvt ++ ;
}
}



void AllpixDigiProcessor::check( LCEvent* ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void AllpixDigiProcessor::end(){

  gsl_rng_free( _rng );
  
  streamlog_out(MESSAGE) << " end()  " << name() 
  << " processed " << _nEvt << " events in " << _nRun << " runs "
  << std::endl ;
  
}
