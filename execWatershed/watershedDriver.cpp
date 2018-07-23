#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if CH_SPACEDIM==3
#define WATERSHED
#else
void THIS_IS_AN_ERROR_MESSAGE(void)
{
  THIS_WILL_ONLY_COMPILES_WHEN_CH_SPACEDIM_IS__3;
}
#endif

#include <petscsnes.h>
#include <iostream>
#include <sys/stat.h>
#include "EBFABView.H"
#include "EBAMRIO.H"
//#include "memtrack.H"
//#Include "memusage.H"
//#include "memtrack.H"
#include "EBDebugDump.H"
#include "DebugOut.H"
#include "ParserBCValue.H"
#include "DebugDump.H"
#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "CH_Attach.H"
#include "EBWatershedUtils.H"
#include "WatershedIBC.H"
#include "EBAMRWatershed.H"
//#include "muParser.h"
#include "EBDebugDump.H"
#include "DebugDump.H"

#include "UsingNamespace.H"


void 
ebamrterrain(const WatershedParameters& a_params,
	     const ProblemDomain& a_coarsestDomain);

int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
    MPI_Init(&a_argc,&a_argv);
#endif
    PetscInitialize(&a_argc,&a_argv, PETSC_NULL, PETSC_NULL);  
    //Scoping trick
    {
	CH_TIMERS("uber_timers");
	CH_TIMER("define_geometry", t1);
	CH_TIMER("run", t2);
	
	//Check for an input file
	char* inFile = NULL;
	if (a_argc > 1)
	{
	    inFile = a_argv[1];
	}
	else
	{
	    pout() << "Usage: <executable name> <inputfile>" << endl;
	    pout() << "No input file specified" << endl;
	    return -1;
	}
	
	//Parse the command line and the input file (if any)
	ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);
	
	ProblemDomain coarsestDomain;
	WatershedParameters params;
    
	//read params from file
	getWatershedParameters(params, coarsestDomain);
    
	//define geometry from given params
	defineWatershedGeometry(params, coarsestDomain);
	
	
	ebamrterrain(params, coarsestDomain);
	
	EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
	ebisPtr->clear();
	
    }//end scoping trick
    PetscFinalize();
#ifdef CH_MPI
    CH_TIMER_REPORT();
    //dumpmemoryatexit();
    MPI_Finalize();
#endif
    return 0;
}

void 
ebamrterrain(const WatershedParameters& a_params,
	     const ProblemDomain& a_coarsestDomain)
{
  
    //CH_TIMERS("ebamrins_driver");
  
  // read inputs
  ParmParse pp;

  Real gconst=1.0;
  pp.query("gravity", gconst);

  // get Richards' equation BCs
  IntVect DomainLoBCType, DomainHiBCType;
  Vector<RefCountedPtr<ParserBCValue> > DomainLoBCValueFunc, DomainHiBCValueFunc;
  {
    Vector<int> itmp(3, 0);
    pp.get("subsurface_domain_Lo_X_type", itmp[0]);
    pp.get("subsurface_domain_Lo_Y_type", itmp[1]);
    pp.get("subsurface_domain_Lo_Z_type", itmp[2]);
    for (int idir=0; idir<3; idir++)
      {
        DomainLoBCType[idir] = itmp[idir];
      }

    pp.get("subsurface_domain_Hi_X_type", itmp[0]);
    pp.get("subsurface_domain_Hi_Y_type", itmp[1]);
    pp.get("subsurface_domain_Hi_Z_type", itmp[2]);
    for (int idir=0; idir<3; idir++)
      {
        DomainHiBCType[idir] = itmp[idir];
      }
    
    Vector<std::string> strtmp(3);
    pp.get("subsurface_domain_Lo_X_func", strtmp[0]);
    pp.get("subsurface_domain_Lo_Y_func", strtmp[1]);
    pp.get("subsurface_domain_Lo_Z_func", strtmp[2]);
    for (int idir=0; idir<3; idir++)
      {
	ParserBCValue* parserBCFuncPtr = new ParserBCValue();
	parserBCFuncPtr->define(RefCountedPtr<ParserFunc>(new ParserFunc(strtmp[idir])),DomainLoBCType[idir]);
	DomainLoBCValueFunc.push_back(RefCountedPtr<ParserBCValue>(parserBCFuncPtr));
      }

    pp.get("subsurface_domain_Hi_X_func", strtmp[0]);
    pp.get("subsurface_domain_Hi_Y_func", strtmp[1]);
    pp.get("subsurface_domain_Hi_Z_func", strtmp[2]);
    for (int idir=0; idir<3; idir++)
      {
	ParserBCValue* parserBCFuncPtr = new ParserBCValue();
	parserBCFuncPtr->define(RefCountedPtr<ParserFunc>(new ParserFunc(strtmp[idir])), DomainHiBCType[idir]);
	//parserBCFuncPtr->setGravityConstant(gconst);
	DomainHiBCValueFunc.push_back(RefCountedPtr<ParserBCValue>(parserBCFuncPtr));
      }
  }
  RefCountedPtr<ParserBCValue>  EBBCValueFunc(NULL);

  if (!a_params.m_includeSurfaceSolver)
  {
      std::string str;
      pp.get("ebbc_sourcevalue_func", str);
      
      ParserBCValue* parserEBBCFuncPtr = new ParserBCValue();
      parserEBBCFuncPtr->define(RefCountedPtr<ParserFunc>(new ParserFunc(str)), 0);
      //parserEBBCFuncPtr->setGravityConstant(gconst);
      EBBCValueFunc = RefCountedPtr<ParserBCValue>(parserEBBCFuncPtr);
  }

  WatershedIBCFactory ibc(DomainLoBCType,
                          DomainHiBCType,
                          DomainLoBCValueFunc,
                          DomainHiBCValueFunc,
			  EBBCValueFunc);


  EBAMRWatershed terrain(a_params, ibc, a_coarsestDomain);

  if (!pp.contains("restart_file"))
  {
      pout() << "starting fresh AMR run" << endl;
      terrain.setupForAMRRun();
  }
  else
  {
      std::string restart_file;
      pp.get("restart_file",restart_file);
      struct stat buffer;   
      if (stat (restart_file.c_str(), &buffer) != 0)
      {
	  MayDay::Error("Restart file does not exist!");
      }
      pout() << " restarting from file " << restart_file << endl;
      terrain.setupForRestart(restart_file);
  }
  
  int maxStep;
  pp.get("max_step", maxStep);

  Real stopTime;
  pp.get("max_time",stopTime);

  Real initialDt;
  pp.get("initial_dt", initialDt);
  terrain.setInitialDt(initialDt);

  terrain.run(stopTime, maxStep);
}
/***************/
