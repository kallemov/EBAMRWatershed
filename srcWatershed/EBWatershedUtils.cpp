#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>
#include "parstream.H"
#include "ParmParse.H"
#include "GeometryShop.H"
#include "EBArith.H"
#include "AllRegularService.H"
#include "PlaneIF.H"
#include "IntersectionIF.H"
#include "UnionIF.H"
#include "TerrainDEMIF.H"
#include "SineChannelIF.H"
#include "EBIndexSpace.H"
#include "EBWatershedUtils.H"

#include "NamespaceHeader.H"


WatershedParameters::WatershedParameters()
{
  m_maxLevel     = 0;
  m_checkpointInterval = -1;
  m_plotInterval       = -1;
  m_maxBoxSize         = 1024;
  m_fillRatio          = 0.7;
  m_blockFactor        = 8;
  m_regridInterval     = -1;
  m_dtGrowFactor       = 1.0;
  
  m_refineThreshold    = 0.1;
  m_verbosity          = 1;
  m_nestingRadius      = 2;
  m_maxDt              = 1.2345678e9;
  
  m_tagShrinkDomain     = 0;
  
  //    m_iterMax             = 100;
  m_orderTimeIntegration = 2;
}
  
WatershedParameters::~WatershedParameters()
{;}


void 
getWatershedParameters(WatershedParameters&   a_params,
		       ProblemDomain&   a_coarsestInputDomain)
{
  // read inputs
  ParmParse ppebamr;
    
  ppebamr.get("init_pressurehead"  ,a_params.m_initial_pressurehead);
  ppebamr.get("solve_overland_flow",a_params.m_includeSurfaceSolver);
  if (a_params.m_includeSurfaceSolver)
  {
      ppebamr.get("ebbc_sourcevalue_func",a_params.m_scalarSurfaceSourceFunc);
  }
  ppebamr.get("include_subsurface_solver",a_params.m_includeSUBSurfaceSolver);
  if (!a_params.m_includeSUBSurfaceSolver && !a_params.m_includeSurfaceSolver)
  {
      MayDay::Error("At least one of the solver must be be included into the solver. Set parameters include_subsurface_solver or/and solve_overland_flow !");
  }

  ppebamr.get("irregular_value_cc",a_params.m_irregular_value_cc);

  Vector<Real> domainLength;
  ppebamr.getarr("domain_length", domainLength,0,SpaceDim);
  for (int idir=0; idir<SpaceDim; idir++)
  {
    a_params.m_domainLength[idir] = domainLength[idir];
  }

  ppebamr.get("max_level"          ,a_params.m_maxLevel);
  ppebamr.get("checkpoint_interval",a_params.m_checkpointInterval);
  ppebamr.get("plot_interval"      ,a_params.m_plotInterval);
  ppebamr.get("max_grid_size"      ,a_params.m_maxBoxSize);
  ppebamr.get("fill_ratio"         ,a_params.m_fillRatio);
  ppebamr.get("block_factor"       ,a_params.m_blockFactor);
  ppebamr.get("regrid_interval"    ,a_params.m_regridInterval);
  ppebamr.query("grow_dt_factor"   ,a_params. m_dtGrowFactor);
  ppebamr.query("max_dt"           ,a_params.m_maxDt);
  ppebamr.query("rerun_smaller_dt" ,a_params.m_useSmallerDt);
  ppebamr.get("surface_output_rate",a_params.m_printOutputRate);
  ppebamr.get("refine_threshold"   ,a_params.m_refineThreshold);
  ppebamr.get("verbosity"          ,a_params.m_verbosity);
  ppebamr.get("nesting_radius"     ,a_params.m_nestingRadius);
  ppebamr.get("tag_buffer"         ,a_params.m_tagBuffer);
  ppebamr.get("refine_all_irregular",a_params.m_refineAllIrreg);
  ppebamr.get("refine_underresolved",a_params.m_refineUnderresolved);
  ppebamr.query("tag_shrinkdomain" ,a_params.m_tagShrinkDomain);
  ppebamr.query("order_time"       ,a_params.m_orderTimeIntegration);
  ppebamr.getarr("ref_ratio"       ,a_params.m_refRatio,0,a_params.m_maxLevel);

  ppebamr.query("tag_function"     ,a_params.m_tagCellsFunc);

  Vector<int> n_cell(SpaceDim);
  ppebamr.getarr("n_cell",n_cell,0,SpaceDim);
  CH_assert(n_cell.size() == SpaceDim);
  IntVect lo = IntVect::Zero;
  IntVect hi;
  for (int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if (n_cell[ivec] <= 0)
        {
          pout() << " bogus number of cells input = " << n_cell[ivec];
          MayDay::Error();
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  // if (ppebamr.contains("solver_tolerance"))
  //   {
  //     ppebamr.get("solver_tolerance", a_params.m_tolerance);
  //   }

  if (ppebamr.contains("domain_periodicity"))
    {
      Vector<int> periodicity(SpaceDim);
      ppebamr.getarr("domain_periodicity", periodicity, 0, SpaceDim);
      bool periodicityArr[SpaceDim];
      for (int idir=0; idir<SpaceDim; idir++)
        {
          periodicityArr[idir] = periodicity[idir]==1;
        }
      a_coarsestInputDomain = ProblemDomain(lo, hi, periodicityArr);
    }
  else
    {
      a_coarsestInputDomain = ProblemDomain(lo, hi);
    }
}


void
defineWatershedGeometry(WatershedParameters&   a_params,
			const ProblemDomain&   a_coarsestDomain,
			int a_whichGeomForced,
			EBIndexSpace* a_ebisPtr)
{
    //there is a bug if we use an anisotropic geometry
    //thus we generate anisotropic grid using a uniform reference grid with unit grid spacing
    //it will require modifications to DEM implicit functions.
    CH_TIME("AMRINSGeometry");

  RealVect coarsestDx;
  RealVect referenceDx = RealVect::Unit;
  for(int idir=0; idir<3;idir++)
  {
      coarsestDx[idir] = a_params.m_domainLength[idir]/Real(a_coarsestDomain.size(idir));
  }

  const int max_level = a_params.m_maxLevel;
  int maxCoarsening = 0;
  
  ProblemDomain finestDomain = a_coarsestDomain;
  for (int ilev = 0; ilev < max_level; ilev++)
  {
      finestDomain.refine(a_params.m_refRatio[ilev]);
      maxCoarsening +=a_params.m_refRatio[ilev]-1;
  }

  RealVect fineDx = referenceDx;
  for (int ilev = 0; ilev < max_level; ilev++)
  {
      fineDx /= a_params.m_refRatio[ilev];
  }

  ParmParse pp;
  RealVect origin;
  Vector<Real> originVect(SpaceDim, 0.0);
  pp.queryarr("domain_origin", originVect, 0, SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
  {
      origin[idir] = originVect[idir]/coarsestDx[idir];
  }


  int whichgeom;
  if (a_whichGeomForced < 0)
    {
      pp.get("which_geom",whichgeom);
    }
  else
    {
      whichgeom = a_whichGeomForced;
    }
  int ebMaxSize = a_params.m_maxBoxSize;

  Real thresholdVoF = 0.0;
  pp.query("threshold_vof",thresholdVoF);
  if (thresholdVoF >= 1) MayDay::Error("thrshd_vof must be less than 1");

  BaseIF *implicitBase=NULL;
  int verbosity = 0;
  if (a_params.m_includeSurfaceSolver)
  {
      
  }
  //3D
  if (whichgeom == 0)
  {
      //allregular
      pout() << "all regular geometry" << endl;
      
      AllRegularService regserv;
      a_ebisPtr->define(finestDomain, origin, fineDx[0], regserv, ebMaxSize, maxCoarsening);
  }
  else if (whichgeom == 1)
  {
      pout() << "hillslope geometry" << endl;
      RealVect slopeNormal;
      vector<Real>  slopeNormalVect(SpaceDim);
      pp.getarr("hillslope_outer_normal",slopeNormalVect, 0, SpaceDim);
      for (int idir = 0; idir < SpaceDim; idir++)
      {
	  //slopeNormal[idir] = slopeNormalVect[idir];
	  slopeNormal[idir] = slopeNormalVect[idir]*coarsestDx[idir];
      }
      
      vector<Real>  slopePointVect(SpaceDim);
      pp.getarr("slope_point",slopePointVect, 0, SpaceDim);
      RealVect slopePoint;
      for (int idir = 0; idir < SpaceDim; idir++)
      {
	  //slopePoint[idir] = slopePointVect[idir];
	  slopePoint[idir] = slopePointVect[idir]/coarsestDx[idir];
      }
      
      bool inside = false;
      
      PlaneIF* slope = new PlaneIF(slopeNormal,slopePoint,inside);
      a_params.m_implicitBaseIF = slope->newImplicitFunction();
      
      implicitBase =slope;
  }
  else if (whichgeom == 2)
  {
      pout() << "vCatchment geometry" << endl;
      
      Real slopeY, slopeX, width;
      pp.get("sides_inclination",slopeY);
      pp.get("channel_inclination",slopeX);
      pp.get("channel_width",width);
      width /= 0.5*coarsestDx[1];
      
      vector<Real>  bottomPointVect(SpaceDim);
      pp.getarr("bottom_center",bottomPointVect, 0, SpaceDim);
      
      RealVect bottomPoint, slopeNormal;
      for (int idir = 0; idir < SpaceDim; idir++)
      {
	  //bottomPoint[idir] = bottomPointVect[idir];
	  bottomPoint[idir] = bottomPointVect[idir]/coarsestDx[idir];
      }
      slopeNormal[0] = slopeX*coarsestDx[0];
      slopeNormal[1] = slopeY*coarsestDx[1];
      slopeNormal[2] = coarsestDx[2];
      bottomPoint[1] -=width;
      
      PlaneIF side1(slopeNormal,bottomPoint, false);
      
      
      slopeNormal[1] *= -1.0;
      bottomPoint[1] +=2.0*width;
      
      PlaneIF side2(slopeNormal,bottomPoint,false);
      
      // create the triangle
      UnionIF triangle(side1,side2);
      
      // create plane to chop tip of triangle
      
      slopeNormal[1] = 0.0;
      bottomPoint[1] -=width;
      
      PlaneIF side3(slopeNormal, bottomPoint,false);
      // create triangle with chopped end
      UnionIF *choppedTriangle =new  UnionIF(triangle,side3);
      
      a_params.m_implicitBaseIF = choppedTriangle->newImplicitFunction();

      implicitBase =choppedTriangle;
  }      
  else if (whichgeom == 3)
  {
      pout() << "DEM geometry" << endl;
      
      //interpType = 1 -> bilinear interpolation
      //interpType = 2 -> bicubic interpolation type1 (this performs cubic interpolation in 1D, matching all value)
      //interpType = 3 -> bicubic interpolation type2 (this performs cubic interpolation in 1D, matching values and derivatives)
      int interpType;
      pp.get("interpType",interpType);
      
      //bottomBuffer is space added below the bathymetry,
      //  (the distance from the deepest spot to the domain box)
      Real bottomElevation;
      pp.get("bottomElevation",bottomElevation);
      
      //verticalScale is used for testing anisotropic vs isotropic geometry, if verticalScale=1.0 then this is the true geometry;
      //   if verticalScale is 100 then all elevations are multiplied by 100...
      Real verticalScale = 1.0;
      pp.query("verticalScale",verticalScale);
      
      //highGround is the elevation given for nodata points with all land neighbors
      //  (useful for higher order interpolation)
      Real highGround = 1.e99;
      pp.query("highGround",highGround);
      Real truncElev = 1.e99;
      
      //   a conversion coeff for elevations values
      Real cellvalue = 1.0;
      pp.query("cellvalue",cellvalue);
      
      //get the file name for the Digital Elevation Model
      std::string demFile;
      pp.get("DEM_file",demFile);
      
      //   a conversion coeff for elevations values
      int edgeType = 1;
      pp.query("edgepointtype", edgeType);

      //we need a true dx here 
      RealVect fineTrueDx = coarsestDx;
      for (int ilev = 0; ilev < max_level; ilev++)
      {
	  fineTrueDx /= a_params.m_refRatio[ilev];
      }
      RealVect trueOrigin = origin;
      trueOrigin *=coarsestDx;
      
      TerrainDEMIF *implicit = new TerrainDEMIF(finestDomain,
						interpType,
						fineTrueDx,
						coarsestDx,
						demFile,
						bottomElevation,
						trueOrigin,
						truncElev,
						highGround,
						verticalScale,
						cellvalue,
						edgeType);
      
      a_params.m_implicitBaseIF = implicit->newImplicitFunction();

      implicitBase =implicit;
  }
  else if (whichgeom == 4)
  {
      pout() << "Sine channel geometry" << endl;
      vector<Real>  PointVect(SpaceDim);
      pp.getarr("channel_reference_point",PointVect, 0, SpaceDim);
      RealVect Point;
      for (int idir = 0; idir < SpaceDim; idir++)
      {
	  Point[idir] = PointVect[idir]/coarsestDx[idir];
      }
      
      Real width;
      pp.get("sinechannel_width",width);
      width /=coarsestDx[1];
      Real slope;
      pp.get("sinechannel_slope",slope);
      slope *=coarsestDx[0]/coarsestDx[2];
      
      Real depth;
      pp.get("sinechannel_depth",depth);
      depth /=coarsestDx[2];
      
      bool inside = false;
      pp.query("sinechannel_inside",inside);
      
      SineChannelIF *channel= new SineChannelIF(Point, width, slope, depth,inside);
      
      a_params.m_implicitBaseIF = channel->newImplicitFunction();
      
      implicitBase =channel;
  }
  else
  {
      //bogus which_geom
      pout() << " bogus which_geom input = "
	     << whichgeom << endl;
      MayDay::Error();
  }
  
  if (!pp.contains("ebis_file"))
  {
      GeometryShop workshop(*implicitBase,verbosity,fineDx, thresholdVoF);
      //this generates the new EBIS
      a_ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, maxCoarsening);
  }
  else
  {
      std::string ebis_file;
      pp.get("ebis_file",ebis_file);
      pout() << " recreating  geometry from file " << ebis_file << endl;
      //define ebis anew from file input
#ifdef CH_USE_HDF5
      HDF5Handle handleIn(ebis_file, HDF5Handle::OPEN_RDONLY);
      a_ebisPtr->define(handleIn);
      handleIn.close();
#endif
  }
  
  delete implicitBase;

  if (pp.contains("ebis_output_file") && !pp.contains("ebis_file"))
  {
      std::string ebis_output_file;
      pp.get("ebis_output_file",ebis_output_file);
#ifdef CH_USE_HDF5      
      HDF5Handle handleOut(ebis_output_file, HDF5Handle::CREATE);
      a_ebisPtr->write(handleOut);
      handleOut.close();
      pout() << "written ebisptr to file: "<< ebis_output_file << endl;
#endif
  }
}
#include "NamespaceFooter.H"
